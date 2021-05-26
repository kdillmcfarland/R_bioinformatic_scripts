"Filter a count table based on minimum number of reads per million in
 minimum number of samples

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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Input parameters:
REQUIRED
  dat = DGEList object of voom normalized data
  min.CPM = numeric minimum counts per million the genes must reach
  name = character string of name for results object. Saved in global environment. Default is
         'dat.filter'
  gene.var = character string of column name containing gene data. Default is 'geneName'. 
             Formerly names
  
REQUIRED (Use one only)
  min.sample = numeric minimum number of samples that must reach min.CPM
  min.pct = numeric minimum percent of samples that must reach min.CPM. In
            whole percent format (10), NOT decimal (0.10)
"

#####

rare.gene.filter <- function(dat, min.CPM, gene.var="geneName", 
                             names="geneName", #depreciated
                             min.sample=NULL, min.pct=NULL,
                             name="dat.filter"){
  ##### Packages #####
  #Check package install
  `%notin%` <- Negate(`%in%`)
  pcks <- c("tidyverse","edgeR")
  pcks.to.install <- pcks[pcks %notin% rownames(installed.packages())]
  if(length(pcks.to.install)){
    print("Please install the following packages.")
    stop(paste0(pcks.to.install)) }
  
  #Load packages
  require(tidyverse)
  require(edgeR)
  
  ##### Check parameters #####
  #Correct input object type?
  if(class(dat) != "DGEList"){ stop("dat object must be a DGEList object") }
  #min sampes or percent only?
  if(!is.null(min.sample) & !is.null(min.pct)){ stop("Please provide only one of min.sample or min.pct") }
  #Replace depreciated names parameter
  if(names != "geneName"){gene.var = names}
  
  ##### Define min number of samples #####
  #Calculate min samples based on percent if provided
  if(!is.null(min.pct)){
    min.sample <- round(nrow(dat$samples)*min.pct/100, digits=0)
  }
  
  ##### List not rare genes #####
  # Convert counts to counts per million
  dat.cpm <- cpm(dat$counts)
  # Calculate number of samples that meet the cutoff per gene
  not.rare.samples <- rowSums(dat.cpm >= min.CPM)
  # List not rare genes to be RETAINED
  not.rare.genes <- names(not.rare.samples[not.rare.samples >= min.sample])

  ##### Filter data to remove rare genes #####
  dat.filter <- dat
  # Filter counts
  dat.filter$counts <- as.data.frame(dat.filter$counts) %>% 
    rownames_to_column() %>% 
    dplyr::filter(rowname %in% not.rare.genes) %>% 
    column_to_rownames()
  
  #If gene info exists, filter as well
  if(!is.null(dat.filter$genes)){
    #Gene name column in gene info table?
    if(gene.var %notin% colnames(dat.filter$genes)){ 
      stop("Gene name varible not present in gene info (dat$genes)") }
      
    # Filter gene key
    dat.filter$genes <- as.data.frame(dat.filter$genes) %>% 
      dplyr::filter(get(gene.var) %in% not.rare.genes)
  }
  
  ##### Save to environment #####
  assign(name, dat.filter, envir = .GlobalEnv)
}