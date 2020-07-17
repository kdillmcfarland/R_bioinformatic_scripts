##### Loop function #####
enrich.fxn <- function(gene.list=NULL,
                       gene.df=NULL, df.group="group",
                       ID.type=c("ENSEMBL","ENTREZ"),
                       category, subcategory=NULL,
                       genome, 
                       basename=NULL, outdir="results/GSEA/"){
  
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
      to.gsea <- gene.df %>% 
        filter(get(df.group) == group.level) %>% 
        dplyr::select(geneName) %>% unlist(use.names = FALSE)
      #Run GSEA and save to results list
      results[[group.level]] <- run.GSEA(to.gsea = to.gsea, 
                                         group.level = group.level,
                                         genome=genome, 
                                         category=category,
                                         subcategory=subcategory)
    }
    
##### Loop through gene lists #####
  } else if(!is.null(gene.list)){
    for(group.level in names(gene.list)){
      print(group.level)
      #Get gene list for each group level
        to.gsea <- gene.list[[group.level]]
        
        #Run GSEA and save to results list 
        results[[group.level]] <- run.GSEA(to.gsea = to.gsea, 
                                           group.level = group.level,
                                           genome=genome, 
                                           category=category,
                                           subcategory=subcategory)
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
    filename <- paste(outdir,"GSEA_",output.name, ".csv", sep="")
  } else if(is.null(basename) & !is.null(subcategory)){
    output.name <- paste(category, gsub(":", ".", subcategory), sep="_")
    filename <- paste(outdir,"GSEA_",output.name, ".csv", sep="")
  } else if(!is.null(basename) & is.null(subcategory)){
    output.name <- paste(basename, category, sep="_") 
    filename <- paste(outdir,"GSEA_",
                      output.name, ".csv", sep="")
  } else{ 
    output.name <- paste(basename, category, gsub(":", ".", subcategory),
                         sep="_") 
    filename <- paste("results/GSEA/GSEA_",
                      output.name, ".csv", sep="")
  }
  
  #Save
  write_csv(results.all, filename)
  
}

##### GSEA function #####
run.GSEA <- function(to.gsea, group.level, ID.type,
                     genome, category, subcategory, ...){
  
  #Convert ENSEMBL IDs if needed
  if(ID.type == "ENSEMBL"){
    #Convert gene list to Entrez ID
    gene.entrez <- clusterProfiler::bitr(to.gsea, fromType="ENSEMBL",
                                         toType=c("ENTREZID","SYMBOL"),
                                         OrgDb=genome)
    
    gene.entrez.list <- gene.entrez$ENTREZID
  } else if(ID.type =="ENTREZ"){
    gene.entrez.list <- to.gsea
  } else{
    stop("Function only allows ENSEMBL or ENTREZ IDs")
  }
  

  #Get database of interest
  if(genome$packageName == "org.Hs.eg.db"){
    
    if(is.null(subcategory)){
      db.species <- as.data.frame(msigdbr(species = "Homo sapiens", 
                                          category = category))
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
             ENSEMBLIDs = paste(ENSEMBL, collapse="/")) %>% 
      ungroup() %>% 
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
      mutate(category=category, subcategory=subcategory,
             #Add columns for group info
             group=group.level, size.group = length(to.gsea)) %>% 
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

