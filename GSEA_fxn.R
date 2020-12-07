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
  or
  gmt_ls = list object with gene ontology data

OPTIONAL
  nperm = number of permutations for P-value calculations. Default is 1000
  name = character to add to output file names. Default NULL
  outdir = Output directory. Default 'results/GSEA/'
  plot = logical if should produce a plot of enrichment sources. Default FALSE
  plot.fdr = Significance cutoff for terms to include in plot. Default 0.05
  plot.groups = List of groups in which to apply FDR cutoff
  plotdir = Output directory for plot. Default the same as outdir except
            in 'figs' instead of 'results'
"

#################

GSEA <- function(gene_list, gmt_file=NULL, gmt_ls=NULL, nperm=1000,
                 name=NULL, outdir="results/GSEA/",
                 plot=FALSE, plot.fdr=0.05, plot.groups=NULL,
                 plotdir='figs/GSEA/'){
  #### Setup ####
  require(tidyverse)
  require(fgsea)
  require(gage)
  
  #Blank list to hold results
  all.results <- list()
  
  #### Data ####
  #Load gene ontology
  if(!is.null(gmt_file)){
    myGO <- fgsea::gmtPathways(gmt_file)
  } else if(!is.null(gmt_ls)){
    myGO <- gmt_ls
  } else {
    stop("Please provide gene ontology data as file for list object.")
  }
  
  #### Loop ####
  #Loop through each list in the gene_list object
  for(genes in names(gene_list)){
    message(genes)
      #Extract 1 gene list
      genes.temp <- gene_list[[genes]]
      #Order by fold change
      genes.temp <- sort(genes.temp, decreasing = TRUE)
      
      #### FGSEA ####
      #Set score type based on fold change
      if(min(genes.temp) < 0 & max(genes.temp) > 0){
        scoreType <- "std"
      } else if(max(genes.temp) <= 0){
        scoreType <- "neg"
      } else if(min(genes.temp) >= 0){
        scoreType <- "pos"
      } else{
        stop("Could not determine score type from fold changes.")
      }
      
      #Run GSEA with fgsea
      fg.result <- fgsea::fgseaSimple(pathways = myGO, 
                                stats = genes.temp,
                                nperm=nperm,
                                #eps=0,
                                scoreType=scoreType) %>% 
        as.data.frame()
      
      #### GAGE ####
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
        dplyr::select(pathway, q.val.x, q.val.y, set.size) %>% 
        pivot_longer(-c(pathway,set.size)) %>% 
        drop_na(value) %>% 
        group_by(pathway, set.size) %>% 
        dplyr::summarise(min.q = min(value, na.rm=TRUE),
                  name = name[value == min.q],
                  .groups = "keep")
      
      #Combine up/down results
      ga.result.format <- ups %>% 
        mutate(gage.FC = "up") %>% 
        filter(pathway %in% filter(gs.signif, name == "q.val.x")$pathway)
      ga.result.format <- downs %>% 
        mutate(gage.FC = "down") %>%
        filter(pathway %in% filter(gs.signif, name == "q.val.y")$pathway) %>% 
        bind_rows(ga.result.format)
      
      #### Combine results ####
      gsea.result <- fg.result %>% 
        dplyr::select(pathway, padj, ES, NES, size) %>% 
        mutate(fgsea.FC = ifelse(NES < 0, "down","up")) %>% 
        dplyr::rename(fgsea.FDR = padj, fgsea.NES = NES, 
               fgsea.ES=ES, fgsea.size=size) %>% 
        full_join(ga.result.format, by="pathway") %>% 
        dplyr::select(pathway:fgsea.FC, q.val, p.geomean,
                      stat.mean, set.size, gage.FC) %>% 
        dplyr::rename(gage.FDR = q.val, gage.geomean = p.geomean, 
                      gage.statmean = stat.mean,
                      gage.size=set.size)  %>% 
        mutate(group = genes)
      
      #### Save ####
      all.results[[genes]] <- gsea.result
  }
  
  #### Format output ####
  #Unlist results into 1 df
  all.results.df <- do.call(rbind.data.frame, all.results)
  rownames(all.results.df) <- NULL
  
  if(!is.null(gmt_file)){
  GO.name <- strsplit(basename(gmt_file), split="[.]")
  obj.name <- paste(GO.name[[1]][[1]], "GSEA.result", sep="_")
  } else {
    obj.name <- "GSEA.result"
  }
  
  assign(obj.name, all.results.df, envir = .GlobalEnv)
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if(!is.null(name)){
    filename <- paste(outdir, obj.name, "_", name,
                      ".csv", sep="")
  } else{
    filename <- paste(outdir, obj.name, ".csv", sep="")
  }
  write_csv(all.results.df, file = filename)
  
  #### Plot ####
  if(plot){
    #Filter to significant results
    if(is.null(plot.groups)){
      to.plot <- all.results.df %>% 
        filter(fgsea.FDR <= plot.fdr)
    } else {
      #First group
      terms.to.plot <- all.results.df %>% 
        filter(fgsea.FDR <= plot.fdr & group %in% plot.groups[[1]]) %>% 
        distinct(pathway) %>% unlist(use.names=FALSE)
      #Remaining groups
      for(i in c(2:length(plot.groups))){
        terms.temp <- all.results.df %>% 
          filter(fgsea.FDR <= plot.fdr & group %in% plot.groups[[i]]) %>% 
          distinct(pathway) %>% unlist(use.names=FALSE)
        
        terms.to.plot <- intersect(terms.to.plot, terms.temp)
      }
      to.plot <- all.results.df %>% 
        filter(pathway %in% terms.to.plot)
    }
  
  
  if(nrow(to.plot > 0)){
  plot.dat <- all.results.df %>% 
    #Significant terms
    filter(pathway %in% to.plot$pathway) %>% 
    #color by significance
    mutate(Significance = ifelse(fgsea.FDR <= 0.01, "FDR < 0.01",
                          ifelse(fgsea.FDR <= 0.05, "FDR < 0.05",
                          ifelse(fgsea.FDR <= 0.1, "FDR < 0.1",
                          ifelse(fgsea.FDR <= 0.2, "FDR < 0.2",
                                 "NS"))))) 
  
  #Enrichment score limits
  plot.lim <- max(abs(plot.dat$fgsea.NES))+0.2
   
  plot <- plot.dat %>%  
    ggplot(aes(reorder(pathway, fgsea.NES), fgsea.NES)) +
    geom_segment(aes(reorder(pathway, fgsea.NES), 
                      xend=pathway, y=0, yend=fgsea.NES)) +
    geom_point(size=3, aes(fill = Significance),
               shape=21, stroke=1) +
    geom_hline(yintercept = 0) +
    
    scale_fill_manual(values=c("FDR < 0.01"="#a50026", 
                               "FDR < 0.05"="#f46d43",
                               "FDR < 0.1"="#fdae61",
                               "FDR < 0.2"="#ffffbf",
                               "NS"="grey")) +
    lims(y=c(-plot.lim,plot.lim)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA", fill = "FGSEA significance") + 
    facet_grid( ~ group) +
    theme_bw()
  
  #### Save plot ####
  dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
  if(!is.null(name)){
    plotname <- paste(plotdir, obj.name, "_", name,
                      ".pdf", sep="")
  } else{
    plotname <- paste(plotdir, obj.name, ".pdf", sep="")
  }
  ggsave(plotname, plot, 
         width = 10+length(unique(to.plot$group))*10, 
         height = length(unique(to.plot$pathway)),
         limitsize = FALSE)
} else{
  message("No significant terms plotted.")
}}
}