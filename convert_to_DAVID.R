"Convert WGCNA gene lists to DAVID format

Used in 'RNAseq_module_fxn.R'

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
"

#################

convert_to_DAVID <- function(voom.dat,
                             mod.net, #From blockwiseModules
                             genes.in.mod,
                             prefix="module",
                             file.out="DAVID.csv"){
  require(tidyverse)
  
  #list modules
  mods <- rownames(voom.dat)
  #calculate maximum module gene list length
  max.mod <- max(table(mod.net$colors))
  #Create data frame with that number of rows
  pDC.david <- data.frame(rowname=1:max.mod)
  
  for (i in 1:length(mods)){
    #Create module name as "module#"
    mod.name = paste(prefix, i, sep=".")
    #Filter gene list to module of interest
    gene.list <- genes.in.mod %>% 
      filter(module == i) %>% 
      dplyr::select(geneName)
    #Calculate the total number of rows that need to be added for nrows to match
    add.genes <- max.mod - nrow(gene.list)
    #Add NAs to fill out gene.list and convert to data frame for merging
    gene.list <- c(gene.list$geneName, rep(NA, times=add.genes))
    gene.list <- as.data.frame(gene.list)
    
    #Combine module lists and rename to mod.name
    pDC.david <- pDC.david %>% 
      bind_cols(gene.list) %>% 
      dplyr::rename(!!quo_name(mod.name) := gene.list)
  }
  
  #Save
  pDC.david %>% 
    dplyr::select(-rowname) %>% 
    write_csv(
      file.out, 
      col_names=TRUE, na="")
}