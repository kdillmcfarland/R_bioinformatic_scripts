"Hypergeometric enrichment of genes in Broad gene sets

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
REQUIRED
ONE OF
  gene.list = List of named vectors as in data.ls['group 1'] = c('gene1','gene2',...) where
              enrichment is assessed separately for each group
  gene.df = data frame with groups in one column and gene IDs in another. Gene IDs must be in geneName
      df.group = Column name in gene.df containing groups to enrich within. Default is 'group'

AND
  category = Character name of Broad gene set to enrich in. One of 'H' or 'C1' through 'C8'
      subcategory = If using a subset of the above gene set. One of
                    'CP' - C2 canonical pathways including BIOCARTA, KEGG, PID, REACTOME
                    'GO' - C5 gene ontology including molecular function (MF), biological process (BP),
                           cellular component (CC)
                    or
                    Colon separated combination of the above with further subsetting as in 
                    'CP:KEGG' or 'GO:BP'
  
  ID.type = Character identifier type for genes in data. One of 'ENTREZ' or 'ENSEBML' or 'SYMBOL'
  genome = Character for genome reference to use. One of 'org.Hs.eg.db', 'org.Mm.eg.db' currently allowed

OPTIONAL
  basename = Character prefix for output file name
  outdir = File path to directory where output is saved. Default is 'results/enrichment/'
   
Example
  enrich.fxn(gene.df=data.df,
             df.group='module',
             category='C5', subcategory='GO',
             ID.type='ENTREZ',
             genome='org.Hs.eg.db', 
             basename='modules', 
             outdir='results/enrichment/')
"

##### Loop function #####
enrich.fxn <- function(gene.list=NULL,
                       gene.df=NULL, df.group="group",
                       species=c("human","mouse"), 
                       category=c("H","C2","C5"), subcategory=NULL,
                       db=NULL, #Add fxn to let user input df with gs_name and gene
                       ID.type=c("symbol","ensembl","entrez")){
  
  ##### Setup #####
  require(clusterProfiler)
  require(msigdbr)
  require(tidyverse)
  
  #Silence warnings
  options(warn=-1)
  
  #Blank holders
  results <- list()
  
  ##### Loop through gene df #####
  if(!is.null(gene.df)){
    #List all possible group levels
    group.list <- unique(gene.df[,df.group])
    
    for(group.level in group.list){
      print(group.level)
      #Get gene list for each group level
      to.enrich <- gene.df %>% 
        filter(get(df.group) == group.level) %>% 
        dplyr::pull(geneName) #Need to make robust to different column names
      #Run enrich and save to results list
      results[[group.level]] <- run.enrich(to.enrich = to.enrich, 
                                           group.level = group.level,
                                           species=species, 
                                           category=category,
                                           subcategory=subcategory,
                                           ID.type=ID.type)
    }
    
    ##### Loop through gene lists #####
  } else if(!is.null(gene.list)){
    for(group.level in names(gene.list)){
      print(group.level)
      #Get gene list for each group level
      to.enrich <- gene.list[[group.level]]
      
      #Run enrich and save to results list 
      results[[group.level]] <- run.enrich(to.enrich = to.enrich, 
                                           group.level = group.level,
                                           species=species, 
                                           category=category,
                                           subcategory=subcategory,
                                           ID.type=ID.type)
    }
    ##### Stop if no genes provided #####
  } else{
    stop("Please provide gene list or data frame.")
  }
  
  ##### Save results #####
  #combine list of df results
  results.all <- dplyr::bind_rows(results, .id = "column_label")
  return(results.all)
}

##### enrich function #####
run.enrich <- function(to.enrich, group.level, 
                       species, category, subcategory, ID.type, ...){
  
  #Get database of interest
  db.species <- as.data.frame(msigdbr(species = species, 
                                      category = category,
                                      subcategory = subcategory))
  
  if(ID.type == "symbol"){
    db.species2 <- dplyr::select(db.species, gs_name, gene_symbol)
  } else if(ID.type == "ensembl"){
    db.species2 <- dplyr::select(db.species, gs_name, ensembl_gene)
  } else if(ID.type == "entrez"){
    db.species2 <- dplyr::select(db.species, gs_name, entrez_gene)
  } else (
    stop("HGNC symbol, ENSEMBL, or ENTREZ ID only.")
  )
  #run enrichment on gene list
  enrich.result <- clusterProfiler::enricher(gene=to.enrich, 
                                             TERM2GENE=db.species2)
  
  #handle no enrichment results
  if (is.null(enrich)){
    enrich.result.clean <- data.frame(
      Description="No enriched terms",
      category=category, 
      group=group.level)
    if (!is.null(subcategory)){
      enrich.result.clean <- enrich.result.clean %>% 
        mutate(subcategory=subcategory)
    }
    return(enrich.result.clean)
    
  }
  else{
    #Format category labels
    db.species.clean <- db.species %>% 
      dplyr::select(gs_cat, gs_subcat, gs_name) %>% 
      dplyr::rename(category=gs_cat, subcategory=gs_subcat, 
                    Description=gs_name) %>% 
      distinct()
    
    #Format results   
    #Format gene column to vector
    enrich.result.clean <- enrich.result@result %>% 
      remove_rownames() %>% 
      arrange(p.adjust, Count) %>% 
      mutate(geneID = strsplit(geneID, split="/")) %>% 
      #Extract values from ratios
      separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
      separate(GeneRatio, into=c("size.overlap.term",
                                 "size.overlap.category"),
               sep="/") %>% 
      mutate_at(vars("size.term","size.category",
                     "size.overlap.term","size.overlap.category"),
                as.numeric) %>% 
      #Calculate k/K
      mutate("k/K"=size.overlap.term/size.term) %>% 
      
      #Add ID columns for database names
      left_join(db.species.clean, by = "Description") %>% 
      #Add columns for group info
      mutate(group=group.level, size.group = length(to.enrich)) %>% 
      #Reorder variables
      dplyr::select(category, subcategory,
                    group, size.group, 
                    size.overlap.category, size.category,
                    Description, size.overlap.term, size.term, `k/K`,
                    p.adjust, qvalue, geneID) %>% 
      arrange(p.adjust)  
    
    return(enrich.result.clean)
  }
}

