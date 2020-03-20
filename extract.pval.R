"Extract p-values from limma results

Creates 2 objects in .GlobalEnv. First is a data frame containing
  Level of variable of interest (generally geneName or module)
  Linear model results including logFC, AveExpr, t, P.Value, adj.P.Val,
  and B
  Variable of interest (group)
  logFC 'up' or 'down' (FC.group)
  
Optional second object is a summary data frame of total genes/modules 
significant at FDR <= 0.05, 0.1, 0.2, 0.3, 0.4, and 0.5

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
                         FC.group = FALSE){
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
      #Calculate total, nonredundant signif genes at different levels
      total.05 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.05) %>% 
        distinct(geneName, FC.group) %>% 
        count(FC.group) %>% 
        rename(`n.05`=n)
      total.1 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.1) %>% 
        distinct(geneName, FC.group) %>% 
        count(FC.group) %>% 
        rename(`n.1`=n)
      total.2 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.2) %>% 
        distinct(geneName, FC.group) %>% 
        count(FC.group) %>% 
        rename(`n.2`=n)
      total.3 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.3) %>% 
        distinct(geneName, FC.group) %>% 
        count(FC.group) %>% 
        rename(`n.3`=n)
      total.4 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.4) %>% 
        distinct(geneName, FC.group) %>% 
        count(FC.group) %>% 
        rename(`n.4`=n)
      total.5 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.5) %>% 
        distinct(geneName, FC.group) %>% 
        count(FC.group) %>% 
        rename(`n.5`=n)
      #Combine results
      gene.tots <- full_join(total.05, total.1, by = "FC.group") %>% 
                   full_join(total.2, by = "FC.group") %>% 
                   full_join(total.3, by = "FC.group") %>% 
                   full_join(total.4, by = "FC.group") %>% 
                   full_join(total.5, by = "FC.group") %>% 
        mutate(group = "total (nonredundant)")
      
      #Summarize signif genes per variable at various levels
      gene.05 <- pval.result %>% 
        filter(adj.P.Val <= 0.05) %>% 
        group_by(group, FC.group, .drop = FALSE) %>% 
        tally() %>% 
        rename(n.05=n)
      
      gene.1 <- pval.result %>% 
        filter(adj.P.Val <= 0.1) %>% 
        group_by(group, FC.group, .drop = FALSE) %>% 
        tally() %>% 
        rename(n.1=n)
      
      gene.2 <- pval.result %>% 
        filter(adj.P.Val <= 0.2) %>% 
        group_by(group, FC.group, .drop = FALSE) %>% 
        tally()  %>% 
        rename(n.2=n)
      
      gene.3 <- pval.result %>% 
        filter(adj.P.Val <= 0.3) %>% 
        group_by(group, FC.group, .drop = FALSE) %>% 
        tally() %>% 
        rename(n.3=n)
      
      gene.4 <- pval.result %>% 
        filter(adj.P.Val <= 0.4) %>% 
        group_by(group, FC.group, .drop = FALSE) %>% 
        tally() %>% 
        rename(n.4=n)
      
      gene.5 <- pval.result %>% 
        filter(adj.P.Val <= 0.5) %>% 
        group_by(group, FC.group, .drop = FALSE) %>% 
        tally()  %>% 
        rename(n.5=n)
      
      #Combine all
      pval.summ <- full_join(gene.05, gene.1, by = c("group", "FC.group")) %>% 
        full_join(gene.2, by = c("group", "FC.group")) %>% 
        full_join(gene.3, by = c("group", "FC.group")) %>% 
        full_join(gene.4, by = c("group", "FC.group")) %>% 
        full_join(gene.5, by = c("group", "FC.group")) %>% 
        ungroup() %>% 
        filter(group != "(Intercept)") %>% 
        bind_rows(gene.tots) %>% 
        mutate(group = fct_relevel(group, "total (nonredundant)", 
                                   after = Inf)) %>% 
        arrange(group)
      
      name2 <- paste(name, "summ", sep=".")
      assign(name2, pval.summ, envir = .GlobalEnv)
    } else if(summary == TRUE){
      #Calculate total, nonredundant signif genes at different levels
      total.05 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.05) %>% 
        distinct(geneName) %>% 
        nrow()
      total.1 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.1) %>% 
        distinct(geneName) %>% 
        nrow()
      total.2 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.2) %>% 
        distinct(geneName) %>% 
        nrow()
      total.3 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.3) %>% 
        distinct(geneName) %>% 
        nrow()
      total.4 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.4) %>% 
        distinct(geneName) %>% 
        nrow()
      total.5 <- pval.result %>% 
        filter(group != '(Intercept)' & adj.P.Val<=0.5) %>% 
        distinct(geneName) %>% 
        nrow()
      #Combine results
      gene.tots <- data.frame(group="total (nonredundant)",
                              n.05=total.05,
                              n.1=total.1,
                              n.2=total.2,
                              n.3=total.3,
                              n.4=total.4,
                              n.5=total.5)
    #Summarize signif genes per variable at various levels
    gene.05 <- pval.result %>% 
      filter(adj.P.Val <= 0.05) %>% 
      group_by(group, .drop = FALSE) %>% 
      tally() %>% 
      rename(n.05=n)
    
    gene.1 <- pval.result %>% 
      filter(adj.P.Val <= 0.1) %>% 
      group_by(group, .drop = FALSE) %>% 
      tally() %>% 
      rename(n.1=n)
    
    gene.2 <- pval.result %>% 
      filter(adj.P.Val <= 0.2) %>% 
      group_by(group, .drop = FALSE) %>% 
      tally()  %>% 
      rename(n.2=n)
    
    gene.3 <- pval.result %>% 
      filter(adj.P.Val <= 0.3) %>% 
      group_by(group, .drop = FALSE) %>% 
      tally() %>% 
      rename(n.3=n)
    
    gene.4 <- pval.result %>% 
      filter(adj.P.Val <= 0.4) %>% 
      group_by(group, .drop = FALSE) %>% 
      tally() %>% 
      rename(n.4=n)
    
    gene.5 <- pval.result %>% 
      filter(adj.P.Val <= 0.5) %>% 
      group_by(group, .drop = FALSE) %>% 
      tally()  %>% 
      rename(n.5=n)
    
    #Combine all
    pval.summ <- full_join(gene.05, gene.1, by = c("group", "FC.group")) %>% 
      full_join(gene.2, by = c("group", "FC.group")) %>% 
      full_join(gene.3, by = c("group", "FC.group")) %>% 
      full_join(gene.4, by = c("group", "FC.group")) %>% 
      full_join(gene.5, by = c("group", "FC.group")) %>% 
      ungroup() %>% 
      filter(group != "(Intercept)") %>% 
      bind_rows(gene.tots) %>% 
      mutate(group = fct_relevel(group, "total (nonredundant)", 
                                 after = Inf)) %>% 
      arrange(group)
    
    name2 <- paste(name, "summ", sep=".")
    assign(name2, pval.summ, envir = .GlobalEnv)
  }
  }