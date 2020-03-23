"Extract p-values from limma results

Creates 2 objects in .GlobalEnv. First is a data frame containing
  Level of variable of interest (generally geneName or module)
  Linear model results including logFC, AveExpr, t, P.Value, adj.P.Val,
  and B
  Variable of interest (group)
  logFC 'up' or 'down' (FC.group)
  
Optional second object is a summary data frame of total genes/modules 
significant at FDR <= list provided

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
  model = Model matrix output from model.matrix( )
  voom.dat = Voom normalized counts data frame, rows as genes and 
             volumns are samples. Output from voom ( ) or 
             voomWithQualityWeights( ) and extracted with 
             as.data.frame(voom.object$E)
  eFit = Empirical Bayes of linear model of interest. Output from
         eBayes(lmFit(voom.dat, model...)) or similar

OPTIONAL
  name = Character string defining name of outputs. Default is 'pval'
  summary = Logical if should calculate summary table with total
            significant genes/modules at several FDR cutoffs. Default 
            is TRUE
  contrasts = Logical if using contrasts matrix. If TRUE, must provide 
              'contrast.mat'. Default is FALSE.
  contrast.mat = Contrast matrix from makeContrasts( ) corresponding to
                 contrasts of variables in the 'model'
  FC.group = Logical if should parse summary for fold change up and down.
            Default is FALSE
  fdr.cutoff = Vector of fdr values to include in summary. Default is
               c(0.05,0.1,0.2,0.3,0.4,0.5)
               
Example
  extract.pval(model = model.interact,
               voom.dat = voom.counts, 
               eFit = efit.data, 
               name = 'pval_mods_2B.interact',
               summary = TRUE,
               contrasts = FALSE,
               FC.group = FALSE)
"

#################
extract.pval <- function(model, voom.dat, eFit, 
                         name="pval",
                         summary=FALSE, 
                         contrast.mat=NULL,
                         contrasts=FALSE,
                         FC.group = FALSE,
                         fdr.cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5)){
  require(tidyverse)
  require(limma)

  #Empty df to hold results
  pval.result <- data.frame()
  
  #List variables of interest
  if(contrasts == TRUE){
    vars <- colnames(contrast.mat)
  } else{
    vars <- colnames(model)
  }
  
  for(var in 1:length(vars)){
    
    pval.temp <- topTable(eFit, coef=var,
                            number=nrow(voom.dat),
                            adjust.method = "BH")
    
    if ("geneName" %in% colnames(pval.temp)){
      pval.temp <- pval.temp %>% 
        mutate(group = vars[var]) %>% 
        # Add categorical var for DE
        mutate(FC.group = ifelse(logFC < 0 , "down",
                                 ifelse(logFC > 0 , "up", NA))) %>% 
        #Convert groups to ordered factors
        mutate(FC.group = factor(FC.group, levels=c("down","up")),
               group = factor(group, levels=vars))
    } else {
      pval.temp <- pval.temp %>% 
        rownames_to_column("geneName") %>% 
        mutate(group = vars[var]) %>% 
        # Add categorical var for DE
        mutate(FC.group = ifelse(logFC < 0 , "down",
                                 ifelse(logFC > 0 , "up", NA))) %>% 
        #Convert groups to ordered factors
        mutate(FC.group = factor(FC.group, levels=c("down","up")),
               group = factor(group, levels=vars))
    }

    pval.result <- bind_rows(pval.result, pval.temp)
  
  }
  assign(name, pval.result, envir = .GlobalEnv)
  
#Add FC group if selected
    if(summary == TRUE & FC.group == TRUE){
      
      total.results <- data.frame()
      group.results <- pval.result %>% 
        count(group, FC.group, .drop = FALSE) %>% 
        select(-n) %>% 
        mutate_if(is.factor, as.character)
      
      for(fdr in fdr.cutoff){
        name.fdr <- paste("fdr",fdr, sep="_")
        #Calculate total, nonredundant signif genes at different levels
        total.temp <- pval.result %>% 
          filter(group != '(Intercept)' & adj.P.Val<=fdr) %>% 
          distinct(geneName, FC.group) %>% 
          count(FC.group, .drop = FALSE) %>% 
          mutate(fdr_group = name.fdr) %>% 
          mutate_if(is.factor, as.character)
        
        total.results <- bind_rows(total.results, total.temp) 
        
        #Summarize signif genes per variable at various levels
        group.temp <- pval.result %>% 
          filter(adj.P.Val <= fdr) %>% 
          count(group, FC.group, .drop = FALSE) %>% 
          select(n) %>% 
          rename(!!name.fdr:="n")
        
        group.results <- bind_cols(group.results, group.temp)
      }
      
      #Combine group and total results
      total.results.pivot <- total.results %>% 
        pivot_wider(names_from = fdr_group, values_from = n) %>% 
        mutate(group = "total (nonredundant)")
      
      result <- group.results %>% 
        filter(group != "(Intercept)") %>% 
        mutate_if(is.factor, as.character) %>% 
        bind_rows(total.results.pivot) 
      
      name2 <- paste(name, "summ", sep=".")
      assign(name2, result, envir = .GlobalEnv)
      
    } else if(summary == TRUE){
          
        total.results <- c()
        group.results <- pval.result %>% 
          count(group, .drop = FALSE) %>% 
          select(-n) %>% 
          mutate_if(is.factor, as.character)
        
        for(fdr in fdr.cutoff){
          name.fdr <- paste("fdr",fdr, sep="_")
          #Calculate total, nonredundant signif genes at different levels
          total.temp <- pval.result %>% 
            filter(group != '(Intercept)' & adj.P.Val<=fdr) %>% 
            distinct(geneName) %>% 
            nrow()
          
          total.results <- c(total.results, total.temp)
          
          #Summarize signif genes per variable at various levels
          group.temp <- pval.result %>% 
            filter(adj.P.Val <= fdr) %>% 
            count(group, .drop = FALSE) %>% 
            select(n) %>% 
            rename(!!name.fdr:="n")
          
          group.results <- bind_cols(group.results, group.temp)
        }
        
        
        #Combine group and total results
        result <- group.results %>% 
          filter(group != "(Intercept)") %>% 
          mutate_if(is.factor, as.character) %>% 
          rbind(c("total (nonredundant)",total.results))
        
        name2 <- paste(name, "summ", sep=".")
        assign(name2, result, envir = .GlobalEnv)
  }
  }