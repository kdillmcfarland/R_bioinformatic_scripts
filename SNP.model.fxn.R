"Run linear models of SNPs across RNA-seq data

Voom normalized gene counts and corresponding SNP data are subset to
samples with both data types input. lmFit linear models of each SNP + any additional model variables and/or sample blocking provided as input vs. expression of a gene of interest. Results saved to disk as 'outfile' value.

#################

Kaitlin Flynn, Kim Dill-Mcfarland, Matt Altman
University of Washington
Copyright (C) 2020
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
  voom.counts = Data frame of voom normalized gene counts with genes as
                rows, samples as columns. Gene names must be in the first
                column, not rownames
  snp.data = Data frame of SNP values for samples in voom.counts. Must be
             0, 1, 2 corresponding to 0/0, 0/1, 1,1. Should also contain
             any additional variables used in snp.model. Variable order
             does not matter
  match.by = Variable name in snp.data that corresponds to voom.counts
             column names. Is used for matching and subsetting these data.
             Default is 'libID'
  gene.no = Character string of ENSEMBL gene number of gene of interst
  snp.model = Character string of additional variables to be added to
              SNP linear model. Must be variable names in snp.data
  
OPTIONAL
  snp.min = Minimum number of samples present at each snp level present in
            the data. Default is 3. If fail, exclusion is defined as 'low
            snp diversity'
  sample.min = Minimum number of samples with SNP call to be used in models. If
               fewer samples than the minimum exist in model matrix, model is not
               run and exclusion is defined as 'low N'. Default is 50
  dupCorr = Logical if models should be blocked. Default is FALSE
  duppCorr.var = Character variable to block samples by. Must be in snp.data. Default is NULL.
  cores = Number of cores to use for parallel computing. Default is 1
  outfile = File path for results file. Default is 'SNP.results.csv'
  
Example
  SNP.model.fxn(voom.counts = as.data.frame(voom.samples$E), 
              snp.data = gene_snps, 
              match.by = 'libID',
              gene.no = 'ENSG00000215182', 
              snp.model = 'gender + age', 
              snp.min = 3, 
              sample.min = 50,
              dupCorr = TRUE,
              dupCorr.var = 'patientID',
              outfile = 'SNP.results.csv', 
              cores = 3)
"

#################

SNP.model.fxn <- function(voom.counts, snp.data, match.by="libID",
                          gene.no, snp.model=NULL, 
                          snp.min=3, sample.min=50,
                          dupCorr=FALSE, dupCorr.var=NULL,
                          outfile="SNP.results.csv", cores=1){
  require(tidyverse)
  require(limma)
  require(data.table)
  require(doParallel)
  require(foreach)
  #Define cores to use
  registerDoParallel(cores=cores)
  
  ##### Subset data #####
  #Keep only libraries with both count and SNP data
  overlap.samp <- snp.data %>% 
    #Force matching name column to be a character, not factor
    mutate_at(match.by, as.character) %>% 
    #Keep samples in SNP data that have count data
    filter(get(match.by) %in% colnames(voom.counts)) %>% 
    #Convert to unnamed vector
    select(match.by) %>% deframe()
  
  #Filter data to data overlap samples
  count.sub <- voom.counts %>% 
    select(1, overlap.samp)
  snp.sub <- snp.data %>%
    select(match.by, everything()) %>% 
    filter(get(match.by) %in% overlap.samp)
  
  ##### Check input data #####
  if(any(colnames(count.sub)[-1] != snp.sub[,1])){
    stop("ERROR: Samples in count and SNP data do not match.")}
  
  ##### Subset to gene of interest ##### 
  count.sub.GOI <- count.sub %>% 
    filter_at(1, all_vars(. == gene.no)) %>% 
    #move gene names to rownames
    rename(rowname=1) %>% 
    column_to_rownames()
  
  snp.sub.GOI <- snp.sub %>% 
    #move matching names to rownames
    column_to_rownames(match.by)
  
  ##### Loop through model for each SNP ##### 
  #List all SNP IDs in data
  snp.list <- colnames(snp.sub.GOI)[grepl("rs[0-9]",
                                          colnames(snp.sub.GOI))]
  #Make empty df to store results
  results <- data.frame()
  
  #LOOP
  results <- rbindlist(fill=TRUE, foreach(i=1:length(snp.list)) %dopar% {
    #Define SNP name
    snp.id <- snp.list[i]
    #List SNP values without NAs
    snp.values <- snp.sub.GOI[!is.na(snp.sub.GOI[,snp.id]),snp.id]

    #Don't run model if SNP does not have multiple levels in data
    if(length(unique(snp.values)) <= 1){
      results.temp <- data.frame(snp=snp.id, 
                                 exclude="no snp diversity")
      results <- bind_rows(results,results.temp)
    } else
      
    #Don't run model if SNP diversity across levels is too low. 
    #Based on snp.min samples per snp level. Exclude = 0 since the previous
      #exclusion takes care of those
    if(any(length(snp.values[snp.values==0]) == c(0:snp.min-1) |
           length(snp.values[snp.values==1]) == c(0:snp.min-1) |
           length(snp.values[snp.values==2]) == c(0:snp.min-1))){
      results.temp <- data.frame(snp=snp.id, 
                                 exclude="low snp diversity")
      results <- bind_rows(results,results.temp)} else
        
    #Don't run model if too few samples with relevant SNP data
    if(length(snp.sub.GOI[!is.na(snp.sub.GOI[,snp.id]),snp.id]) < sample.min){
      results.temp <- data.frame(snp=snp.id, 
                                 exclude="low N")
      results <- bind_rows(results,results.temp)
    } else{
      #Else run model
      
      # Define model matrix
      if(!is.null(snp.model)){
        curM <- paste("~snp.sub.GOI[,snp.id]", snp.model, sep="+")
        curMM <- model.matrix(as.formula(curM), data=snp.sub.GOI)
      } else{
        curMM <- model.matrix(~snp.sub.GOI[,snp.id], data=snp.sub.GOI)
      }
      
      #Models WITHOUT blocking
      #If all samples have SNP data and remain in model,
      #Run model on data as is
      if (dupCorr==FALSE & nrow(curMM) == nrow(snp.sub.GOI))
      {
        fitAgs <- lmFit(count.sub.GOI, curMM)
        fitAgs <- eBayes(fitAgs)
        results.temp <- topTable(fitAgs, coef=2,
                                 number=nrow(count.sub.GOI),
                                 sort.by="none",
                                 adjust.method="none") %>% 
          mutate(snp=snp.id)
      }
      
      #If 1+ samples are missing SNP data,
      #Remove these samples and run model on subset data
      if (dupCorr==FALSE & nrow(curMM) < nrow(snp.sub.GOI))
      {
        #List samples with missing SNP data
        remIndex <- which(!(rownames(snp.sub.GOI) %in%
                              rownames(curMM)))
        #Remove samples and run model
        fitAgs <- lmFit(count.sub.GOI[,-remIndex], curMM)
        fitAgs <- eBayes(fitAgs)
        results.temp <- topTable(fitAgs, coef=2,
                                 number=nrow(count.sub.GOI),
                                 sort.by="none",
                                 adjust.method="none") %>% 
          mutate(snp=snp.id)
      }
      
      #Models WITH blocking
      #If all samples have SNP data and remain in model,
      #Run model on data as is
      if (dupCorr==TRUE & nrow(curMM) == nrow(snp.sub.GOI))
      {
        dupCorr.value <- duplicateCorrelation(count.sub.GOI, curMM, 
                                      block=snp.sub.GOI[,dupCorr.var])
        fitAgs <- lmFit(count.sub.GOI, curMM,
                        block=snp.sub.GOI[,dupCorr.var],
                        correlation=dupCorr.value$consensus.correlation)
        fitAgs <- eBayes(fitAgs)
        results.temp <- topTable(fitAgs, coef=2,
                                 number=nrow(count.sub.GOI),
                                 sort.by="none",
                                 adjust.method="none") %>% 
          mutate(snp=snp.id)
      }
      
      #If 1+ samples are missing SNP data,
      #Remove these samples and run model on subset data
      if (dupCorr==TRUE & nrow(curMM) < nrow(snp.sub.GOI))
      {
        #List samples with missing SNP data
        remIndex <- which(!(rownames(snp.sub.GOI) %in%
                              rownames(curMM)))
        #Remove samples and run model
        dupCorr.value <- duplicateCorrelation(count.sub.GOI[,-remIndex],
                                            curMM, 
                                            block=snp.sub.GOI[-remIndex,
                                                              dupCorr.var])
        fitAgs <- lmFit(count.sub.GOI[,-remIndex], curMM,
                        block=snp.sub.GOI[-remIndex,dupCorr.var],
                        correlation=dupCorr.value$consensus.correlation)
        fitAgs <- eBayes(fitAgs)
        results.temp <- topTable(fitAgs, coef=2,
                                 number=nrow(count.sub.GOI),
                                 sort.by="none",
                                 adjust.method="none") %>% 
          mutate(snp=snp.id)
      }
        
      }
      })
  #Save final results to disk
  write_csv(results, outfile)
}