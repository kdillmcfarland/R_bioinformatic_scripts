##### Loop function #####
enrich.fxn <- function(gene.list=NULL,
                       gene.df=NULL, df.group="group",
                       category, subcategory=NULL,
                       ID.type=NULL,
                       genome, 
                       basename=NULL, outdir="results/enrichment/"){
  
##### Setup #####
  require(clusterProfiler)
  require(msigdbr)
  require(tidyverse)
  require(plyr)
  if(genome$packageName == "org.Hs.eg.db"){
    require(org.Hs.eg.db)}
  if(genome$packageName == "org.Mm.eg.db"){
    require(org.Mm.eg.db)}
  
  #Silence warnings
  options(warn=-1)

  #Blank holders
  results <- list()
  
##### Loop through gene df #####
  if(!is.null(gene.df)){
    #List all possible group levels
    group.list <- unique(gene.df[,df.group]) %>% 
      dplyr::select(all_of(df.group)) %>% unlist(use.names = FALSE)
    
    for(group.level in group.list){
      print(group.level)
      #Get gene list for each group level
      to.enrich <- gene.df %>% 
        filter(get(df.group) == group.level) %>% 
        dplyr::select(geneName) %>% unlist(use.names = FALSE)
      #Run enrich and save to results list
      results[[group.level]] <- run.enrich(to.enrich = to.enrich, 
                                         group.level = group.level,
                                         genome=genome, 
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
                                           genome=genome, 
                                           category=category,
                                           subcategory=subcategory,
                                           ID.type=ID.type)
    }
##### Stop if no genes provided #####
  } else{
    stop("Please provide gene list or data frame.")
  }
  
##### Save results #####
  dir.create(outdir, showWarnings = FALSE)
  
  #combine list of df results
  results.all <- plyr::ldply (results, data.frame) %>% 
    dplyr::select(-'.id')
  
  #Make filename
  if(is.null(basename) & is.null(subcategory)){ 
    output.name <- category 
    filename <- paste(outdir,"enrich_",output.name, ".csv", sep="")
  } else if(is.null(basename) & !is.null(subcategory)){
    output.name <- paste(category, gsub(":", ".", subcategory), sep="_")
    filename <- paste(outdir,"enrich_",output.name, ".csv", sep="")
  } else if(!is.null(basename) & is.null(subcategory)){
    output.name <- paste(basename, category, sep="_") 
    filename <- paste(outdir,"enrich_",
                      output.name, ".csv", sep="")
  } else{ 
    output.name <- paste(basename, category, gsub(":", ".", subcategory),
                         sep="_") 
    filename <- paste(outdir, "enrich_",
                      output.name, ".csv", sep="")
  }
  
  #Save
  write_csv(results.all, filename)
  
}

##### enrich function #####
run.enrich <- function(to.enrich, group.level, 
                     genome, category, subcategory, ID.type, ...){
  
  #Convert ENSEMBL IDs if needed
  if(ID.type == "ENSEMBL"){
    #Convert gene list to Entrez ID
    gene.entrez <- clusterProfiler::bitr(to.enrich, fromType="ENSEMBL",
                                         toType=c("ENTREZID","SYMBOL"),
                                         OrgDb=genome)
    
    gene.entrez.list <- gene.entrez$ENTREZID
  } else if(ID.type =="ENTREZ"){
    gene.entrez <- clusterProfiler::bitr(to.enrich, fromType="ENTREZID",
                                         toType=c("ENSEMBL","SYMBOL"),
                                         OrgDb=genome)
    gene.entrez.list <- to.enrich
  } else{
    stop("Function only allows ENSEMBL or ENTREZ IDs")
  }
  

  #Get database of interest
  if(genome$packageName == "org.Hs.eg.db"){
    
    #No subcategory
    if(is.null(subcategory)){
      db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                          category = category))
    } else
    # Combine all CP subs
    if(subcategory == "CP"){
      db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                          category = "C2",
                                          subcategory = "CP:BIOCARTA")) %>% 
        bind_rows(as.data.frame(msigdbr(species = "Homo sapiens", 
                                        category = "C2",
                                        subcategory = "CP:KEGG"))) %>% 
        bind_rows(as.data.frame(msigdbr(species = "Homo sapiens", 
                                        category = "C2",
                                        subcategory = "CP:PID"))) %>% 
        bind_rows(as.data.frame(msigdbr(species = "Homo sapiens", 
                                        category = "C2",
                                        subcategory = "CP:REACTOME")))
    } else
      # Combine all GO subs
      if(subcategory == "GO"){
        db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                            category = "C5",
                                            subcategory = "GO:MF")) %>% 
          bind_rows(as.data.frame(msigdbr(species = "Homo sapiens", 
                                          category = "C5",
                                          subcategory = "GO:BP"))) %>% 
          bind_rows(as.data.frame(msigdbr(species = "Homo sapiens", 
                                          category = "C5",
                                          subcategory = "GO:CC")))
      } else if(subcategory=="GO:MF"){
        db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                            category = "C5",
                                            subcategory = "GO:MF"))
    } else if(subcategory=="GO:BP"){
      db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                          category = "C5",
                                          subcategory = "GO:BP"))
    } else if(subcategory=="GO:CC"){
      db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                          category = "C5",
                                          subcategory = "GO:CC")) %>% 
        bind_rows(as.data.frame(msigdbr(species = "Homo sapiens", 
                                        category = "C5",
                                        subcategory = "GO:BP"))) %>% 
        bind_rows(as.data.frame(msigdbr(species = "Homo sapiens", 
                                        category = "C5",
                                        subcategory = "GO:CC")))
    } else {
      db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                          category = category,
                                          subcategory = subcategory))}
    
  } else if(genome$packageName == "org.Mm.eg.db"){
    
    if(is.null(subcategory)){
      db.species <- as.data.frame(msigdbr(species = "Mus musculus",
                                          category = category)) 
    } else{
      db.species <- as.data.frame(msigdbr(species = "Mus musculus",
                                          category = category,
                                          subcategory = subcategory)) 
    }
  } else{
    stop("Function only available for human and mouse genomes.")
  }
  
  #run enrichment on gene list
  enrich <- enricher(gene=gene.entrez.list, 
                     TERM2GENE=dplyr::select(db.species, gs_name,
                                             entrez_gene))
  
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
    #Extract results
    enrich.result <- enrich@result %>% 
      remove_rownames() %>% 
      arrange(p.adjust, Count)
    
    #Create group names for entrez+number genes ID'ed
    ## Use to separate list of entrez IDs if > 1 exist for a given term
    pivot_names <- c()
    for (i in 1:max(enrich.result$Count)){
      pivot_names[[i]] <- paste("entrez", i, sep="")
    }
    
    #Format category labels
    db.species.clean <- db.species %>% 
      dplyr::select(gs_cat, gs_subcat, gs_name) %>% 
      dplyr::rename(category=gs_cat, subcategory=gs_subcat, 
                    Description=gs_name) %>% 
      distinct()
    #Format results   
    enrich.result.clean <- enrich.result %>% 
      #Separate entrez ID lists
      separate(geneID, into=as.character(pivot_names), sep="/") %>% 
      pivot_longer(all_of(as.character(pivot_names)), names_to = "rep", 
                   values_to = "ENTREZID") %>% 
      drop_na(ENTREZID) %>% 
      #Match entrez IDs to gene IDs
      left_join(gene.entrez, by="ENTREZID") %>% 
      
      #Combine lists into single columns, sep by /
      group_by_at(vars(ID:Count)) %>% 
      dplyr::summarize(ENTREZIDs = paste(ENTREZID, collapse="/"),
             SYMBOLs = paste(SYMBOL, collapse="/"),
             ENSEMBLIDs = paste(ENSEMBL, collapse="/"),
             .groups="drop") %>% 
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
                    p.adjust, qvalue, ENTREZIDs:ENSEMBLIDs) %>% 
      arrange(p.adjust)  
    
    return(enrich.result.clean)
  }
}

