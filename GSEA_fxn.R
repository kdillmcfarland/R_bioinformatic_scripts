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
  fdr = numeric FDR cutoff for significance
  
"

#################

GSEA <- function(gene_list, gmt_file, fdr=0.05){
  require(tidyverse)
  require(fgsea)
  require(gage)
  
  #Blank list to hold results
  all.results <- list()
  #Load gene ontology
  myGO <- fgsea::gmtPathways(gmt_file)
  
  #Loop through each list in the gene_list object
  for(genes in names(gene_list)){
      #Extact 1 gene list
      genes.temp <- gene_list[[genes]]
      #Order by fold change
      genes.temp <- sort(genes.temp, decreasing = TRUE)
      
      #Run GSEA with fgsea
      fg.result <- fgsea::fgsea(pathways = myGO, 
                                stats = genes.temp) %>% 
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
      gs.signif <- full_join(ups, downs, by="pathway") %>% 
        select(pathway, q.val.x, q.val.y) %>% 
        pivot_longer(-pathway) %>% 
        group_by(pathway) %>% 
        summarise(min.q = min(value),
                  name = name[value == min.q])
      
      ga.result.format <- ups %>% 
        filter(pathway %in% filter(gs.signif, name == "q.val.x")$pathway) %>% 
        mutate(gage.FC = "up")
      ga.result.format <- downs %>% 
        filter(pathway %in% filter(gs.signif, name == "q.val.y")$pathway) %>% 
        bind_rows(ga.result.format) %>% 
        mutate(gage.FC = ifelse(is.na(gage.FC), "down", gage.FC))
      
      #Combine two methods.
      gsea.result <- fg.result %>% 
        select(pathway, padj, NES) %>% 
        mutate(fgsea.FC = ifelse(NES < 0, "down","up")) %>% 
        rename(fgsea.FDR = padj, fgsea.NES = NES) %>% 
        full_join(ga.result.format) %>% 
        select(pathway:fgsea.FC, q.val, stat.mean, gage.FC) %>% 
        rename(gage.FDR = q.val, gage.NES = stat.mean) %>% 
        mutate(concordant = ifelse(fgsea.FDR <= 0.05 & gage.FDR <= 0.05 & 
                            fgsea.FC == gage.FC, "signif.concordant", NA)) %>% 
        mutate(group = genes)
      
      #save to object
      all.results[[genes]] <- gsea.result
  }
  
  #Unlist results into 1 df
  all.results.df <- do.call(rbind.data.frame, all.results)
  rownames(all.results.df) <- NULL
  assign("GSEA.result", all.results.df, envir = .GlobalEnv)
}