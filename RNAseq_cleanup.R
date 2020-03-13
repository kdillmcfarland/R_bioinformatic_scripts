"Clean default BRI RNA-seq data outputs

Includes 'Final Annotation.xlsx', 'combined_summary-data.csv', and
'combined_counts.csv'

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
  data.dir = Filepath to directory containg all files
                         
OPTIONAL
  anno.cols = Columns to extract from annotation file
              Default is c('Proj-Sub...1', 'Donor ID', 
              'Study Group', 'Date Collected',
              'Sample Id', 'cDNA sampleId','library sampleId')
  summ.cols = Columns to extract from summary file
              Default is c('libId', 'fastq_total_reads', 'total_reads',
              'total_sequences', 'total_counts', 'median_cv_coverage',
              'mapped_reads_w_dups', 'predicted_sex')
"

#################

clean.RNAseq <- function(
  data.dir=NULL,
  #Default lists of column names
  anno.cols=c("Proj-Sub...1", "Donor ID", 
              "Study Group", "Date Collected",
              "Sample Id", "cDNA sampleId","library sampleId"),
  summ.cols=c("libId", "fastq_total_reads", "total_reads",
              "total_sequences", "total_counts", "median_cv_coverage",
              "mapped_reads_w_dups", "predicted_sex")
                     ){
  require(tidyverse)
  require(readxl)
  require(stringr)
  require(lubridate)
  
  ##### Final annotation #####
  #Find file path
  final.anno.file <- dir(data.dir, pattern = "*Final Annotation.xlsx", 
                    full.names = TRUE) %>% 
    gsub("//", "/", .)
  #Load file
  final.anno <- read_excel(final.anno.file)
  
  #Clean
  final.anno.format <- final.anno %>% 
    #Select column of interest
    select(anno.cols) %>% 
    #Rename ID columns
    rename_if(startsWith(names(.), "Proj"), ~"projID") %>% 
    rename_if(startsWith(names(.), "Donor"), ~"donorID") %>% 
    rename("sampID" = "Sample Id",
           "cDNAID" = "cDNA sampleId") %>% 
    rename_if(startsWith(names(.), "library"), ~"libID") %>% 
    #Rename metadata columns
    rename("group" = "Study Group") %>% 
    rename_if(startsWith(names(.), "Date"), ~"date") %>% 
    #format date data
    mutate(date = as.Date(date))
  
  ##### Combined summary #####
  #Find file path
  summ.file <- dir(data.dir, pattern = "*combined_summary-data.csv", 
                         full.names = TRUE) %>% 
    gsub("//", "/", .)
  #Load file
  summ <- read_csv(summ.file)
  
  #Extract data of interest
  summ.format <- summ %>% 
    select(summ.cols) %>% 
    separate(libId, into = c("libID","flow.cell"), sep="_")

  ##### Combine metadata #####
  meta.clean <- full_join(final.anno.format, summ.format, by="libID") %>% 
    arrange(donorID, libID) 
  #Save to environment
  assign("meta.clean", meta.clean, envir = .GlobalEnv)
  
  ##### Counts table #####
  #Find file path
  count.file <- dir(data.dir, pattern = "*combined_counts.csv", 
                   full.names = TRUE) %>% 
    gsub("//", "/", .)
  #Load file
  counts <- read_csv(count.file)
  
  #Remove flowcell ID from column names
  counts.format <- counts %>% 
    set_names(~gsub("_.*", "", names(counts))) %>% 
    select(geneName, meta.clean$libID)
  #Save to environment
  assign("counts", counts.format, envir = .GlobalEnv)
}
            