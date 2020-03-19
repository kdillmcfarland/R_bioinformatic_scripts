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
  min.CPM = minimum counts per million the genes must reach
  name = name of results object. Saved in global environment. Default is
         'dat.filter'
  
REQUIRED (Use one only)
  min.sample = Numeric minimum number of samples that must reach min.CPM
  min.pct = Numeric minimum percent of samples that must reach min.CPM. In
            whole percent format (10), not decimal (0.10)
"

#####

rare.gene.filter <- function(dat, min.CPM,
                             min.sample=NULL, min.pct=NULL,
                             name="dat.filter"){
  require(tidyverse)
  require(edgeR)
  
  if(class(dat) != "DGEList"){
    stop("dat object must be a DGEList")
  }
  
  ##### Define min number of samples #####
  if(!is.null(min.pct)){
    min.sample <- round(nrow(dat$samples)*min.pct/100, digits=0)
  }
  
  ##### List not rare genes #####
  # Convert counts to counts per million
  dat.cpm <- cpm(dat$counts)
  # Calculate number of samples that meet the cutoff per gene
  not.rare.samples <- rowSums(dat.cpm >= min.CPM)
  # List not rare genes to be RETAINED
  not.rare.genes <- dat[not.rare.samples >= min.sample, ]$genes$geneName
  
  ##### Filter data to remove rare genes #####
  # Filter counts
  dat.filter <- dat
  dat.filter$counts <- as.data.frame(dat.filter$counts) %>% 
    rownames_to_column("geneName") %>% 
    dplyr::filter(geneName %in% not.rare.genes) %>% 
    column_to_rownames("geneName")
  # Filter gene key
  dat.filter$genes <- as.data.frame(dat.filter$genes) %>% 
    dplyr::filter(geneName %in% not.rare.genes)
  
  ##### Save to environment #####
  assign(name, dat.filter, envir = .GlobalEnv)
}