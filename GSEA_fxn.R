"Run GSEA with fgsea and gage

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
  gene_list = named list object with named numeric vectors of gene symbols and logFC
  gmt_file = character giving full path the gene ontology GMT file

"

#################

GSEA <- function(gene_list, gmt_file, name=NULL){
  require(tidyverse)
  require(fgsea)
  require(gage)
  
  #Blank list to hold results
  all.results <- list()
  #Load gene ontology
  myGO <- fgsea::gmtPathways(gmt_file)
  
  #Loop through each list in the gene_list object
  for(genes in names(gene_list)){
    message(genes)
      #Extact 1 gene list
      genes.temp <- gene_list[[genes]]
      #Order by fold change
      genes.temp <- sort(genes.temp, decreasing = TRUE)
      
      #Run GSEA with fgsea
      fg.result <- fgsea::fgsea(pathways = myGO, 
                                stats = genes.temp,
                                eps=0) %>% 
        as.data.frame()
      
      #Run GSEA with gage
      ga.result <- gage::gage(genes.temp, gsets=myGO)
      
      #Extract FC up
      ups = as.data.frame(ga.result$greater) %>% 
        rownames_to_column("pathway")
      #Extract FC down
      downs = as.data.frame(ga.result$less) %>% 
        rownames_to_column("pathway")
      #Combine most significant q value for each term
      gs.signif <- full_join(ups, downs, by=c("pathway", "set.size")) %>% 
        select(pathway, q.val.x, q.val.y, set.size) %>% 
        pivot_longer(-c(pathway,set.size)) %>% 
        group_by(pathway, set.size) %>% 
        summarise(min.q = min(value),
                  name = name[value == min.q],
                  .groups = "keep")
      
      ga.result.format <- ups %>% 
        filter(pathway %in% filter(gs.signif, name == "q.val.x")$pathway) %>% 
        mutate(gage.FC = "up")
      ga.result.format <- downs %>% 
        filter(pathway %in% filter(gs.signif, name == "q.val.y")$pathway) %>% 
        bind_rows(ga.result.format) %>% 
        mutate(gage.FC = ifelse(is.na(gage.FC), "down", gage.FC))
      
      #Combine two methods.
      gsea.result <- fg.result %>% 
        select(pathway, padj, ES, NES, size) %>% 
        mutate(fgsea.FC = ifelse(NES < 0, "down","up")) %>% 
        rename(fgsea.FDR = padj, fgsea.NES = NES, 
               fgsea.ES=ES, fgsea.size=size) %>% 
        full_join(ga.result.format, by="pathway") %>% 
        select(pathway:fgsea.FC, q.val, p.geomean, stat.mean, set.size, gage.FC) %>% 
        rename(gage.FDR = q.val, gage.geomean = p.geomean, gage.statmean = stat.mean,
               gage.size=set.size)  %>% 
        mutate(group = genes)
      
      #save to object
      all.results[[genes]] <- gsea.result
  }
  
  #Unlist results into 1 df
  all.results.df <- do.call(rbind.data.frame, all.results)
  rownames(all.results.df) <- NULL
  
  GO.name <- strsplit(basename(gmt_file), split="[.]")
  obj.name <- paste(GO.name[[1]][[1]], "GSEA.result", sep="_")
  assign(obj.name, all.results.df, envir = .GlobalEnv)
  
  dir.create("results/GSEA_FoldChange/", showWarnings = FALSE)
  if(!is.null(name)){
    filename <- paste("results/GSEA_FoldChange/", obj.name, "_", name,
                      ".csv", sep="")
  } else{
    filename <- paste("results/GSEA_FoldChange/", obj.name, ".csv", sep="")
  }
  write_csv(all.results.df, path = filename)
}