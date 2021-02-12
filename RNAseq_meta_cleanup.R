RNAseq.meta.clean <- function(data.dir=NULL, 
                              out.file="data_clean/cleaning.metrics.csv",
                              trim=TRUE,
                              align=TRUE,
                              bam=TRUE,
                              bam.filter=TRUE,
                              count=TRUE
                              ){
  require(tidyverse)

  summ.all <- data.frame(libID=NA)
  
  #### Raw and adapter trimmed ####
  if(trim){
    # Raw and adapter trimmmed seqs
    trim.files <- list.files(data.dir, pattern="*settings", 
                             all.files=FALSE, full.names=TRUE) %>% 
      gsub("//", "/", .)
    
    trim.summ <- data.frame()
    
    for (file in trim.files){
      libID.temp <- basename(file) %>% 
        gsub(".settings", "", .)
      
      trim.summ.temp <- read_table(file, col_names = FALSE,
                                   col_types = cols()) %>% 
        filter(startsWith(X1, 'Total number of read pairs') |
                 startsWith(X1, 'Number of retained reads')) %>% 
        separate(X1, into=c("metric","value"), sep=": ") %>% 
        mutate(libID = libID.temp) %>% 
        pivot_wider(names_from = "metric", values_from = "value") %>% 
        rename(raw="Total number of read pairs",
               trim="Number of retained reads") %>% 
        mutate_at(vars(raw, trim), as.numeric) %>% 
        mutate(raw=raw*2)
      
      trim.summ <- bind_rows(trim.summ, trim.summ.temp)
    }
    
    # Combine into results df
    summ.all <- full_join(summ.all, trim.summ, by = "libID")
  }
  
  #### Alignment ####
  if(align){
    align.file <- list.files(data.dir, pattern="summary.alignment", 
                             all.files=FALSE, full.names=TRUE) %>% 
      gsub("//", "/", .)
    
    align.lib <- read_tsv(align.file, col_names = FALSE,
                          col_types = cols()) %>% 
      filter(grepl("bam",X1)) %>% 
      distinct(X1)
      
    # Align summary
    align.summ <- read_tsv(align.file, col_names = FALSE,
                           col_types = cols()) %>% 
      #Get sample name from filename
      mutate(libID = factor(X1, levels=align.lib$X1),
             libID = basename(as.character(libID)),
             libID = gsub("_Aligned.sortedByCoord.out.bam", "", libID)) %>% 
      fill(libID) %>% 
      #Separate data to column
      separate(X1, into=c("h","i"), sep=" \\+ 0 ", fill="right") %>% 
      drop_na(i) %>% 
      #Recode data types (f)
      separate(i, into=c("i"), sep="[(]", extra="drop") %>% 
      mutate(i = fct_recode(factor(i),
                            to.be.aligned="in total ",
                            secondary.align="secondary",
                            chimeric.align="supplementary",
                            PCR.dups="duplicates",
                            align="mapped ",
                            paired="paired in sequencing",
                            R1.paired="read1", R2.paired="read2",
                            align.paired="properly paired ",
                            both.align.paired= "with itself and mate mapped",
                            one.align.paired="singletons " ,
                            both.align.paired.diffCHR="with mate mapped to a different chr",
                            both.align.paired.diffCHR.mapq="with mate mapped to a different chr ")) %>% 
      pivot_wider(names_from = "i", values_from = "h") %>% 
      mutate_at(vars(-libID), as.numeric)
    
    # Combine into results df
    summ.all <- full_join(summ.all, align.summ, by = "libID")
  } 
  
  #### BAM summary ####
  if(bam){
    bam1 <- list.files(data.dir, pattern="bam.metrics", 
                       all.files=FALSE, full.names=TRUE) %>% 
      gsub("//", "/", .)
    
    bam.colnames <- read_delim(bam1, delim = " ", col_names=FALSE)[7,1] %>% 
      separate(X1, into=as.character(c(1:30)), sep="\t") %>% 
      unlist(use.names = FALSE)
    
    bam.summ <- read_delim(bam1, delim = " ", col_names=FALSE, comment = "#",
                           col_types = cols()) %>%
      #Get sample name from file name
      mutate(libID = factor(X1, levels=align.lib$X1),
             libID = basename(as.character(libID)),
             libID = gsub("_Aligned.sortedByCoord.out.bam", "", libID)) %>% 
      fill(libID) %>% 
      #Sep data columns
      separate(X1, into = bam.colnames, sep="\t", fill="right") %>% 
      #Remove histogram columns
      drop_na(CODING_BASES) %>% 
      #Keep only data rows
      filter(!startsWith(as.character(PF_BASES), 'PF')) %>% 
      #Make numeric
      mutate_at(vars(-libID), as.numeric) %>% 
      #Calculate perc align
      mutate(PCT_PF_ALIGNED = PF_ALIGNED_BASES/PF_BASES) %>% 
      select(libID, PCT_PF_ALIGNED, everything())
    
    # Combine into results df
    summ.all <- full_join(summ.all, bam.summ, by = "libID")
  } 
  
  #### Filtered paired alignment ####
  if(bam.filter){
    align.file2 <- list.files(data.dir, pattern="summary.align.filter", 
                             all.files=FALSE, full.names=TRUE) %>% 
      gsub("//", "/", .)
    
    align.lib2 <- read_tsv(align.file2, col_names = FALSE,
                           col_types = cols()) %>% 
      filter(grepl("bam",X1)) %>% 
      distinct(X1)
    
    # Align summary
    align.summ2 <- read_tsv(align.file2, col_names = FALSE,
                            col_types = cols()) %>% 
      #Get sample name from filename
      mutate(libID = factor(X1, levels=align.lib2$X1),
             libID = basename(as.character(libID)),
             libID = gsub("_filter_paired.bam", "", libID)) %>% 
      fill(libID) %>% 
      #Separate data to column
      separate(X1, into=c("h","i"), sep=" \\+ 0 ", fill="right") %>% 
      drop_na(i) %>% 
      #Recode data types (f)
      separate(i, into=c("i"), sep="[(]", extra="drop") %>% 
      mutate(i = fct_recode(factor(i),
                            to.be.aligned="in total ",
                            secondary.align="secondary",
                            chimeric.align="supplementary",
                            PCR.dups="duplicates",
                            align="mapped ",
                            paired="paired in sequencing",
                            R1.paired="read1", R2.paired="read2",
                            align.paired="properly paired ",
                            both.align.paired= "with itself and mate mapped",
                            one.align.paired="singletons " ,
                            both.align.paired.diffCHR="with mate mapped to a different chr",
                            both.align.paired.diffCHR.mapq="with mate mapped to a different chr ")) %>% 
      pivot_wider(names_from = "i", values_from = "h") %>% 
      mutate_at(vars(-libID), as.numeric)
    
    # Combine into results df
    summ.all <- full_join(summ.all, align.summ2, by = "libID")
  }

  #### Gene counts ####
  if(count){
    count.file <- list.files(data.dir, pattern="featurecounts", 
                             all.files=FALSE, full.names=TRUE) %>% 
      gsub("//", "/", .)
    
    count.summ <- read_tsv(count.file, col_types = cols()) %>% 
      pivot_longer(-Status) %>% 
      #Get sample names from file
      mutate(libID = basename(name),
             libID = gsub("_filter.paired.bam", "", libID)) %>% 
      select(-name) %>% 
      #Keep rows with non-zero values
      filter(value !=0) %>% 
      #Transpose
      pivot_wider(names_from = Status, values_from = value)
    
    # Combine into results df
    summ.all <- full_join(summ.all, count.summ, by = "libID")
  }
  
  #### Save ####
  #remove columns that are all blank or 0
  summ.all.filter <- summ.all %>%
    drop_na(libID) %>% 
    mutate_at(vars(-libID), as.numeric) %>% 
    select_if(~sum(!is.na(.)) > 0) %>% 
    select_if(~sum(.!= 0) > 0)
  
  write_csv(summ.all.filter, file=out.file)
}