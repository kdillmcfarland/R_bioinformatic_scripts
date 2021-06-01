"Run linear models with or without random effects and pairwise kinship

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
DATA
REQUIRED 
  dat = EList object containing normalized gene expression (E), sample metadata (targets), and gene
        information (genes). Output by voom() or voomWithQualityWeights()
  OR
  counts = data frame or matrix with log2 CPM. Rows are genes, columns are libraries
  meta = data frame with library metadata. Must contain 'libID' with values that match 
         columns in counts
  gene.info = data frame with gene annotations. Must contain 'geneName' with values that 
              match rows in counts
  
MAIN MODEL
REQUIRED
  x.var = character vector of x-variables to use in model
  ptID = character string of variable name of IDs to match expression, meta, and kinship data. 
         Default is 'FULLIDNO'
  OR
  full.model = character string of model to use such as '~ variable'. Do not include random effects
               as this is added from ptID based on lmekin, lme, or lm usage
  ptID = character string of variable name of IDs to match expression, meta, and kinship data. 
         Default is 'FULLIDNO'
  
OPTIONAL
  kin = if running models with kinship, a numeric matrix with pairwise kinship values. Rows must
        be in the same order as dat and have rownames
  co.var = character vector of co-variates to use in model
  interaction = logical if should include interaction between first two x.var. Default is FALSE
         
OTHER MODEL INFO
OPTIONAL
  lm, lme = logical if should run corresponding simple linear model or linear mixed effects model 
            without kinship for comparison to full model. Defaults are FALSE
  lme.pairwise = logical if should run pairwise comparisons within multiple levels of lme
                 Default is FALSE
  subset.var = character string of variable name to run subsets of data. Must be 1 variable. Useful
               for pairise contrast comparisons
  subset.lvl = character string of level of subset.var to subset to
      For example, subset.var = 'condition', subset.lvl = 'MEDIA'
  subset.genes = character vector of genes to test. If not given, function runs all genes in dat

OTHER, OPTIONAL
  p.method = Method of FDR adjustment. Default is 'BH'
  outdir = character string of output directory. Default is 'results/gene_level/', 
  name = character string of prefix for output file names. Default is 'lme.results'
  processors = Numeric for parallel processors. Default is 1
  
"

#################

lmekin.loop <- function(dat=NULL, counts=NULL, meta=NULL, gene.info=NULL,
                        kin=NULL, x.var, ptID="FULLIDNO",
                        co.var=NULL, interaction=FALSE, 
                        lm=FALSE, lme=FALSE, lme.pairwise=FALSE,
                        full.model = NULL,
                        subset.var = NULL, subset.lvl = NULL, subset.genes = NULL,
                        outdir="results/gene_level/", name="lme.results",
                        processors=1, p.method="BH"){
#Log start time
old <- Sys.time()

##### Packages #####
#Check package install
`%notin%` <- Negate(`%in%`)
pcks <- c("tidyverse","data.table","limma","lme4","coxme","emmeans",
          "broom","car","foreach","doParallel")
pcks.to.install <- pcks[pcks %notin% rownames(installed.packages())]
if(length(pcks.to.install)){
  print("Please install the following packages.")
  stop(paste0(pcks.to.install)) }

#Load packages
##Table manipulation
require(tidyverse, quietly = TRUE,warn.conflicts = FALSE)
require(data.table, quietly = TRUE,warn.conflicts = FALSE)
##Work with input EList object
require(limma, quietly = TRUE,warn.conflicts = FALSE)
##Linear mixed effects models
library(lme4, quietly = TRUE,warn.conflicts = FALSE)
library(emmeans, quietly = TRUE,warn.conflicts = FALSE)
##Linear mixed effects models with kinship
require(coxme, quietly = TRUE,warn.conflicts = FALSE)
##Extract model results
library(broom, quietly = TRUE,warn.conflicts = FALSE)
library(car, quietly = TRUE,warn.conflicts = FALSE)
##Parallel for loops
require(foreach, quietly = TRUE,warn.conflicts = FALSE)
require(doParallel, quietly = TRUE,warn.conflicts = FALSE)

###### Parallel ###### 
#setup parallel processors
registerDoParallel(processors)

###### Check common input parameter errors #####
if(is.null(subset.var) & !is.null(subset.lvl)){
  stop("Sample subsetting has been selected. Please also provide subset.var")}
if(!is.null(subset.var) & is.null(subset.lvl)){
  stop("Sample subsetting has been selected. Please also provide subset.lvl")}
if(!is.null(full.model)){
  if(grepl("[|]", full.model)){
    stop("full.model should not include random effects such as (1|ptID). Please correct")} }

###### Data #####
print("Load data")

#If data are NOT a voom EList, create a mock version
if(is.null(dat)) {
  dat.format <- list()
  
  #Expression data
  ##Move rownames to column if exist
  ##Order columns as in metadata and genes as in gene.info
  if(rownames(counts)[1]!=1){
    counts.format <- as.data.frame(counts) %>% 
      rownames_to_column() %>% 
      select(rowname, all_of(meta$libID)) %>% 
      arrange(match(rowname, gene.info$geneName)) %>% 
      column_to_rownames()
  } else {
    counts.format <- as.data.frame(counts) %>% 
      rename_if(is.character, ~"rowname")%>% 
      select(rowname, all_of(meta$libID)) %>% 
      arrange(match(rowname, gene.info$geneName)) %>% 
      column_to_rownames()
  }
  
  #Metadata
  ##Remove samples not in expression data
  meta.format <- meta %>% 
    filter(libID %in% colnames(counts.format))
  
  #Put in list
  dat.format$E <- counts.format
  dat.format$targets <- meta
  dat.format$genes <- gene.info
} else {
  dat.format <- dat
}

#Format data
#If has rownames, move into df
if(is.numeric(dat.format$E[,1])){
  dat.format$E <- as.data.frame(dat.format$E) %>% 
    rownames_to_column()
} else {
#Rename 1st column
  colnames(dat.format$E)[1] <- "rowname"
}

###### Subset to variable of interest if selected ######
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
  #Combine expression data (E) and sample metadata (targets)
  to.model <- dat.subset$E %>% 
    pivot_longer(-rowname, names_to = "libID", values_to = "expression") %>% 
    inner_join(dat.subset$targets, by="libID") %>% 
    #Remove samples missing kinship
    filter(get(ptID) %in% colnames(kin))
  
  #Compute number of samples to run in models
  rna.no <- dat.subset$targets %>% 
    distinct(get(ptID)) %>% nrow()
  kin.no <- to.model %>% 
    distinct(get(ptID)) %>% nrow()
  
  message(paste(rna.no-kin.no, "individuals missing kinship data. Running models on", 
                kin.no))
}else{
  #Combine expression data (E) and sample metadata (targets)
  to.model <- dat.subset$E %>% 
    pivot_longer(-rowname, names_to = "libID", values_to = "expression") %>% 
    inner_join(dat.subset$targets, by="libID")
  
  #Compute number of samples to run in models
  rna.no <- to.model %>% 
    distinct(get(ptID)) %>% nrow()
  
  message(paste("No kinship provided. Running models on",  rna.no, "individuals"))
}


###### Run models ######
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
  
  #### Simple LM models, if selected #####
  #Run linear model without kinship
  #Place holder LM results
  p.lm <- NaN
  sigma.lm <- 0
  results.lm <- NULL

  if(lm){
    #Make LM formula. as.formula does not work 
    if(!is.null(full.model)) {
      model.lm <- paste("expression", full.model, sep="")
    } else if(interaction){
      model.lm <- paste("expression ~ ", paste(x.var, collapse=" * "), " + ", 
                     paste(co.var, collapse=" + "), 
                     sep="")
    } else {
      model.lm <- paste("expression ~ ", paste(x.var, collapse=" + "), " + ", 
                     paste(co.var, collapse=" + "), 
                     sep="")
    }
    
    #Wrap model run in error catch to allow loop to continue even if a single model fails
    tryCatch({
      #Fit model
      fit.lm <- lm(model.lm, data=to.model.gene)
          p.lm <- tidy(fit.lm)
          sigma.lm <- sigma(fit.lm)
        
      #Extract results 
      results.lm <- data.frame(
        model = rep("lm", nrow(p.lm)),    #Label model as lm
        gene = rep(gene, nrow(p.lm)),     #gene name
        variable = p.lm$term,             #variables in model
        pval = p.lm$p.value,              #P-value
        sigma = rep(sigma.lm, nrow(p.lm)))#sigma
    }, error=function(e){})
  }
  
  #### Simple LME models, if selected #####
  #Make LME formula. as.formula does not work 
  if(!is.null(full.model)) {
    model <- paste("expression", full.model, " + ", 
                   "(1|",ptID,")",
                   sep="") 
    } else if(interaction){
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
  
  #Place holder LME results
  p.lme <- NaN
  sigma.lme <- 0
  results.lme <- NULL
  
  if(lme){
    tryCatch({
      #Fit LME model
      fit.lme <- lmer(model, data=to.model.gene)
          #Estimate P-value
          p.lme <- tidy(Anova(fit.lme))
          #Calculate sigma
          sigma.lme <- sigma(fit.lme)
          
      #Extract results
      results.lme <- data.frame(
        model = rep("lme", nrow(p.lme)),    #Label model as lme
        gene = rep(gene, nrow(p.lme)),      #gene name
        variable = p.lme$term,              #variables in model
        pval = p.lme$p.value,               #P-value
        sigma = rep(sigma.lme, nrow(p.lme)))#sigma
    }, error=function(e){})
  }
  
  ##### Kinship model ######
  #Place holder LMEKIN results
  p.kin <- NaN
  sigma.kin <- 0
  results.kin <- NULL
  
  if(!is.null(kin)){
  tryCatch({
    #Fit LMEKIN model
        fit.kin <- lmekin(as.formula(model), data=to.model.gene, varlist=as.matrix(kin))
            #Calulate stats
            beta <- fit.kin$coefficients$fixed
            nvar <- length(beta)
            nfrail <- nrow(fit.kin$var) - nvar
            se <- sqrt(diag(fit.kin$var)[nfrail + 1:nvar])
            t <- beta/se
            p.kin <- signif(1 - pchisq((t)^2, 1), 2)
            sigma.kin <- fit.kin$sigma
        
        #Extract results
        results.kin <- data.frame(
            model = rep("lmekin", length(p.kin)),  #Label model as lmekin
            gene = rep(gene, length(p.kin)),       #gene name
            variable = names(p.kin),               #variables in model
            pval = p.kin,                          #P-value
            sigma = rep(sigma.kin, length(p.kin))) #sigma
        }, error=function(e){})
  }
  
  #### Pairwise within lme ####
  p.pair <- NaN
  sigma.pair <- 0
  results.pair <- NULL
  
  if(lme.pairwise){
    #For each x.var
    for(x.i in x.var){
      fit.lsmeans <- lsmeans(fit.lme, as.formula(paste("pairwise~",x.i,sep=""))) 
      
      results.pair <- as.data.frame(fit.lsmeans$contrasts) %>% 
        select(contrast, p.value) %>% 
        rename(variable=contrast, pval=p.value) %>% 
        mutate(model="lsmeans", gene=gene)
    }
  }
  
  #### Combine results #####
  #All models for this gene
  results <- results.lm %>% 
    bind_rows(results.lme) %>% 
    bind_rows(results.kin) %>% 
    bind_rows(results.pair)
  
  #This gene to all previous gene results
  fit.results <- rbind(results, fit.results) 
  })

#### Calculate FDR ####
fit.results.fdr <- fit.results %>% 
  #Within model and variable
  group_by(model, variable) %>% 
  mutate(FDR=ifelse(model != "lsmeans",p.adjust(pval, method=p.method),
                    NA)) %>% 
  ungroup() %>% 
  #Add identifier name to allow for easy combination with other lmekin.fxn() outputs
  mutate(group=name) %>% 
  dplyr::select(group, everything())

#### Save ####
print("Saving results")
#Create output directory if not present
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
#Define filename
filename <- paste(outdir, name, ".model.results.csv.gz", sep="")
#Save to csv. gz compress to save space
write.table(fit.results.fdr, sep=",", row.names=FALSE, col.names=TRUE,
            file=gzfile(filename))

###### Fin ###### 
print("All models complete")
#Print total time to run
new <- Sys.time() - old
print(new)
}
