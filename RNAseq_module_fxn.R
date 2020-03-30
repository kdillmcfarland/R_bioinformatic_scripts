"Create WGCNA modules and format outputs

#################

Kim Dill-Mcfarland
University of Washington, kadm@uw.edu
Copyright (C) 2020 Kim Dill-Mcfarland
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
  voom.dat = EList object output by voom( ) or voomWithQualityWeights( )
  genes.signif = Character vector of genes in include in modules
  minModuleSize = Minimum module size. Default is 20
  deepSplit = Level of deep split in module building tree. Default is 3.
              0-4 allowed
              
REQUIRED (one of)
  Rsq.min = Minimum R-squared allowed to select soft thresholding in
            pickSoftThreshold( ).
  sft.value = Override min R-squared and set soft thresholding to a value
  
OPTIONAL
  nThread = Number of parallel processors to use. Default is 1
  basename = Character to prepend to output names. Default is 'basename'
  outdir = Subdirectory name within results/ and figs/ in which to
           save ouputs

Example
  make.modules(voom.dat = dat.voom,
             genes.signif = genes.signif,
             Rsq.min = 0.895,
             minModuleSize = 20,
             deepSplit = 3,
             nThread = 4,
             basename = 'P259.2_allGenes',
             outdir = NULL)
"

#################
make.modules <- function(voom.dat,
                         genes.signif,
                         Rsq.min = NULL,
                         sft.value = NULL,
                         minModuleSize = 20,
                         deepSplit = 3,
                         nThread=1, basename="basename",
                         outdir=NULL){
  require(tidyverse)
  
  ##### Load data ##### 
  voom.dat <- voom.dat
  
  ##### Subset significant genes #####
  voom.signif <- voom.dat
  voom.signif$E <- as.data.frame(voom.signif$E) %>% 
    rownames_to_column() %>% 
    filter(rowname %in% genes.signif) %>% 
    column_to_rownames()
  
  voom.signif$genes <- voom.signif$genes %>% 
    rownames_to_column() %>% 
    filter(rowname %in% genes.signif) %>% 
    column_to_rownames()
  
  ##### Soft-thresholding #####
  allowWGCNAThreads(nThreads=nThread)
  #Calculate
  sft <- pickSoftThreshold(t(voom.signif$E), 
                           powerVector=c(1:30), verbose=5,
                           networkType = "signed")
  #Select threshold 
  if(!is.null(Rsq.min)){
    sft.select <- as.data.frame(sft$fitIndices) %>% 
      filter(SFT.R.sq  >= Rsq.min)
    power.t <- min(sft.select$Power)
  } else if(!is.null(sft.value)){
    sft.select <- as.data.frame(sft$fitIndices) %>% 
      filter(Power == sft.value)
    power.t <- unique(sft.select$Power)
  } else{
    stop("Please select min R-squared or soft threshold.")
  }
  
  #Plot
  sft.plot <- as.data.frame(sft$fitIndices) %>% 
    mutate(fit = -sign(slope)*SFT.R.sq) %>% 
    pivot_longer(c(fit,mean.k.), names_to = "name", values_to = "value") %>% 
    
    ggplot(aes(x=Power, y=value, label=Power)) +
    geom_text(size=3) +
    facet_wrap(~name, scales = "free_y",
               labeller=labeller(name=
               c(fit="Scale free topology model fit,\nsigned R^2",
                 mean.k.="Mean connectivity"))) +
    theme_classic() +
    labs(y="", x="Soft threshold (power)") +
    geom_hline(data=data.frame(name=c("fit","mean.k."), 
                               cutoff=c(sft.select$SFT.R.sq[1],
                                        sft.select$mean.k.[1])),
               aes(yintercept = cutoff), color="red") +
    geom_vline(xintercept = power.t, color="red") 
  
  ##### Build modules #####
  mod.net <- blockwiseModules(t(voom.signif$E),
                              power=power.t, 
                              networkType="signed",
                              TOMType="signed",
                              maxBlockSize=500,
                              minModuleSize=minModuleSize,
                              deepSplit=deepSplit, 
                              numericLabels=TRUE,
                              saveTOMFileBase="TOM-blockwise",
                              nthreads=nThread)
  
  mods <- as.data.frame(mod.net$colors) %>% 
    rownames_to_column("geneName") %>% 
    dplyr::rename(module = "mod.net$colors") %>% 
    left_join(voom.signif$genes, by="geneName") %>% 
    #add leading 0 to module names for correct sorting of factor
    mutate(module.char = ifelse(module <= 9, 
                                paste("0", module, sep=""),
                                module)) %>% 
    #Add color var
    mutate(mod.color = labels2colors(mod.net$colors))
  
  ##### Mean module expression #####
  voom.mods <- mods %>% 
    #Combine count and module data
    dplyr::select(geneName, module.char) %>% 
    left_join(rownames_to_column(voom.signif$E, "geneName"),
              by="geneName") %>% 
    
    #Calculate mean by module
    group_by(module.char) %>% 
    summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
    
    #Add basename to module names
    mutate(module.char=paste("module", basename, module.char, sep="_")) %>% 
    column_to_rownames("module.char")
  
  ##### Save results to environ #####
  assign("sft.select", sft.select[1,], envir = .GlobalEnv)
  
  ##### Write results #####
  #Create directories
  dir.basename <- paste("module_", basename, 
                        "_deepSplit", deepSplit, 
                        "_minMod", minModuleSize, sep="")
  if(!is.null(outdir)){
    dir.results <- paste("results/", outdir, "/", 
                         dir.basename, sep="")
    dir.figs <- paste("figs/", outdir, "/", 
                      dir.basename, sep="")
  } else{
    dir.results <- paste("results/", dir.basename, sep="")
    dir.figs <- paste("figs/", dir.basename, sep="")
  }
  
  dir.create(dir.results, showWarnings = FALSE)
  dir.create(dir.figs, showWarnings = FALSE)
  
  #Sft threshold plot
  ggsave(plot=sft.plot, width=10, height=4,
         filename = paste(dir.figs, "/SFT_thresholding_power", power.t, 
                          ".png", sep=""))
  
  #Voom counts
  write_csv(rownames_to_column(voom.mods, "module"), col_names = TRUE,
            paste(dir.results, "/", basename,
                  "_mod_voom_counts.csv", sep=""))

  #Genes in module
  write_csv(mods, col_names = TRUE,
            paste(dir.results, "/",basename,
                  "_genes_in_mod.csv", sep=""))
  
  #Genes in module, DAVID format
  source("https://raw.githubusercontent.com/kdillmcfarland/R_bioinformatic_scripts/master/convert_to_DAVID.R")
  convert_to_DAVID(voom.dat=voom.mods,
                   mod.net=mod.net,
                   genes.in.mod=mods,
                   prefix=basename,
                   file.out=paste("results/", dir.basename, "/",basename,
                                  "_genes_in_mod_DAVID.csv", sep=""))
}
