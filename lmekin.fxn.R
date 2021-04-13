"Run linear mixed effects models

#################

Kim Dill-Mcfarland
University of Washington, kadm@uw.edu
Copyright (C) 2021 Kim Dill-Mcfarland
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Input parameters:
REQUIRED
  dat = EList object containing normalized gene expression (E), sample metadata (targets), and gene
        information (genes). Output by voom or voomWithQualityWeights
  x.var = character vector of x-variables to use in model
  ptID = character string of variable name of IDs to match expression, meta, and kinship data. 
         Default is 'FULLIDNO'
  p.method = Method of FDR adjustment. Default is 'BH'
  
OPTIONAL
  kin = if running models with kinship, a numeric matrix with pairwise kinship values. Must be in the
        same order as dat
  co.var = character vector of co-variates to use in model
  interaction = logical if should include interaction between first two x.var. Default is FALSE
  lm, lme = logical if should run corresponding simple linear model or linear mixed effects model 
            without kinship for comparison to full model. Defaults are FALSE
  subset.var = character string of variable name to run subsets of data. Must be 1 variable. Useful
               for contrasts model
  subset.lvl = character string of level of subset.var to subset to
      For example, subset.var = 'condition', subset.lvl = 'MEDIA'
  subset.genes = character vector of genes to test. If not given, function runs all genes in dat
  outdir = character string of output directory. Default is 'results/gene_level/', 
  name = character string of prefix for output file names Default is 'lme.results'
  processors = Numeric for parallel processors. Default is 1
  
"

#################

lmekin.loop <- function(dat, kin=NULL, x.var, ptID="FULLIDNO",
                        co.var=NULL, interaction=FALSE, 
                        lm=FALSE, lme=FALSE,
                        subset.var = NULL, subset.lvl = NULL, subset.genes = NULL,
                        outdir="results/gene_level/", name="lme.results",
                        processors=1, p.method="BH"){
old <- Sys.time()
require(tidyverse, quietly = TRUE,warn.conflicts = FALSE)
require(data.table, quietly = TRUE,warn.conflicts = FALSE)
require(limma, quietly = TRUE,warn.conflicts = FALSE)
library(broom, quietly = TRUE,warn.conflicts = FALSE)
library(lme4, quietly = TRUE,warn.conflicts = FALSE)
library(car, quietly = TRUE,warn.conflicts = FALSE)
require(coxme, quietly = TRUE,warn.conflicts = FALSE)
require(foreach, quietly = TRUE,warn.conflicts = FALSE)
require(doParallel, quietly = TRUE,warn.conflicts = FALSE)

###### Parallel ###### 
#setup parallel processors
registerDoParallel(processors)

###### Data #####
print("Load data")

dat.format <- dat
#If rownames, move into df
if(is.numeric(dat$E[,1])){
  dat.format$E <- as.data.frame(dat.format$E) %>% 
    rownames_to_column()
} else {
#Rename 1st column
  colnames(dat.format$E)[1] <- "rowname"
}

###### Subset if selected ######
dat.subset <- dat.format

#Subset samples
if(!is.null(subset.var)){
  dat.subset$targets <- dat.subset$targets %>% 
    filter(get(subset.var) == subset.lvl)
  
  dat.subset$E <- as.data.frame(dat.subset$E) %>% 
    dplyr::select(rowname, all_of(dat.subset$targets$libID))
}

#Subset genes
if(!is.null(subset.genes)){
  dat.subset$E <- as.data.frame(dat.subset$E) %>% 
    filter(rowname %in% subset.genes)
}

###### Format data for modeling ####
if(!is.null(kin)){
  to.model <- dat.subset$E %>% 
    pivot_longer(-rowname, names_to = "libID", values_to = "expression") %>% 
    inner_join(dat.subset$targets, by="libID") %>% 
    #Remove samples missing kinship
    filter(get(ptID) %in% colnames(kin))
  
  rna.no <- dat.subset$targets %>% 
    distinct(get(ptID)) %>% nrow()
  kin.no <- to.model %>% 
    distinct(get(ptID)) %>% nrow()
  
  message(paste(rna.no-kin.no, "individuals missing kinship data. Running models on", 
                kin.no))
}else{
  to.model <- dat.subset$E %>% 
    pivot_longer(-rowname, names_to = "libID", values_to = "expression") %>% 
    inner_join(dat.subset$targets, by="libID")
  
  rna.no <- to.model %>% 
    distinct(get(ptID)) %>% nrow()
  
  message(paste("No kinship provided. Running models on",  rna.no, "individuals"))
}


###### Run model with kinship ######
print("Run models")

#create blank df to hold results
fit.results <- data.frame()

#Loop through each gene
fit.results <- rbindlist(fill=TRUE, foreach(i=1:nrow(dat.subset$E)) %dopar% {
  #### Prepare data ####
  #Get gene name
  gene <- dat.subset$E[i,1]
message(gene)
  
  #Filter data to gene
  to.model.gene <- to.model %>% 
    filter(rowname == gene) %>% 
    arrange(ptID)
  
  #Make model
  if(interaction){
    model <- paste("expression ~ ", paste(x.var, collapse=" * "), " + ", 
                   paste(co.var, collapse=" + "), " + ", 
                   "(1|",ptID,")",
                   sep="")
  } else {
    model <- paste("expression ~ ", paste(x.var, collapse=" + "), " + ", 
                   paste(co.var, collapse=" + "), " + ", 
                   "(1|",ptID,")",
                   sep="")
  }
  
  #### Basic models, if selected #####
  p.lm <- NaN
  sigma.lm <- 0
  results.lm <- NULL
  
  if(lm){
    if(interaction){
      model.lm <- paste("expression ~ ", paste(x.var, collapse=" * "), " + ", 
                     paste(co.var, collapse=" + "), 
                     sep="")
    } else {
      model.lm <- paste("expression ~ ", paste(x.var, collapse=" + "), " + ", 
                     paste(co.var, collapse=" + "), 
                     sep="")
    }
    
    tryCatch({
      fit.lm <- lm(model.lm, data=to.model.gene)
          p.lm <- tidy(fit.lm)
          sigma.lm <- sigma(fit.lm)
          
      results.lm <- data.frame(
        model = rep("lm", nrow(p.lm)),
        gene = rep(gene, nrow(p.lm)),
        variable = p.lm$term, 
        pval = p.lm$p.value,
        sigma = rep(sigma.lm, nrow(p.lm)))
    }, error=function(e){})
  }
  
  p.lme <- NaN
  sigma.lme <- 0
  results.lme <- NULL
  
  if(lme){
    tryCatch({
      fit.lme <- lmer(model, data=to.model.gene)
          p.lme <- tidy(Anova(fit.lme))
          sigma.lme <- sigma(fit.lme)
          
      results.lme <- data.frame(
        model = rep("lme", nrow(p.lme)),
        gene = rep(gene, nrow(p.lme)),
        variable = p.lme$term, 
        pval = p.lme$p.value,
        sigma = rep(sigma.lme, nrow(p.lme)))
    }, error=function(e){})
  }
  
  ##### Kinship model ######
  p.kin <- NaN
  sigma.kin <- 0
  results.kin <- NULL
  
  if(!is.null(kin)){
  tryCatch({
        fit.kin <- lmekin(as.formula(model), data=to.model.gene, varlist=as.matrix(kin))
            beta <- fit.kin$coefficients$fixed
            nvar <- length(beta)
            nfrail <- nrow(fit.kin$var) - nvar
            se <- sqrt(diag(fit.kin$var)[nfrail + 1:nvar])
            t <- beta/se
            p.kin <- signif(1 - pchisq((t)^2, 1), 2)
            sigma.kin <- fit.kin$sigma
            
        results.kin <- data.frame(
            model = rep("lmekin", length(p.kin)),
            gene = rep(gene, length(p.kin)),
            variable = names(p.kin), 
            pval = p.kin,
            sigma = rep(sigma.kin, length(p.kin)))
        }, error=function(e){})
  }
      
  #### Combine results #####
  results <- results.lm %>% 
    bind_rows(results.lme) %>% 
    bind_rows(results.kin) 
      
  fit.results <- rbind(results, fit.results) 
  })

#### Calculate FDR ####
fit.results.fdr <- fit.results %>% 
  group_by(model, variable) %>% 
  mutate(FDR=p.adjust(pval, method=p.method)) %>% 
  ungroup() %>% 
  mutate(group=name) %>% 
  dplyr::select(group, everything())

#### Save ####
print("Saving results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
filename <- paste(outdir, name, ".model.results.csv.gz", sep="")
write.table(fit.results.fdr, sep=",", row.names=FALSE, col.names=TRUE,
            file=gzfile(filename))

###### Fin ###### 
print("All models complete")
new <- Sys.time() - old
print(new)
}