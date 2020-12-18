"Create STRING network of for gene list

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
  genes = Vector of genes to plot in STRING. Can be HGNC symbol, ENSEMBL ID, 
          or ENTREZ ID
  version = Character for STRING database version. Default is '11'
  score_threshold = Numeric for minimum combined score to include. Default is
                    400 which is considered medium confidence
  
OPTIONAL
Coloring by enrichment terms
  enrichment = Matrix/data frame output by hypergeo_enricher function. Necessary
               columns are Description, size.overlap.term, p.adjust, and 
  ID = Character for ID type of gene list. One of c('SYMBOLs','ENTREZIDs','ENSEMBLIDs')
  size.overlap.term = Gene overlap minimum for enrichment terms to include in plot.
                      Default is 2 
  p.adjust = FDR maximum for enrichment terms to include in plot. Default is 0.2
  colors = Vector of colors for manual scales. Default it to use ggplot colors.
Saving plot
  basename = Character to prepend to output names.
  outdir = Character for output directory for plot
  height, width = Plot dimensions in inches
"

#################
string.plot <- function(genes, version="11", score_threshold=400,
                        enrichment=NULL, size.overlap.term=2, p.adjust=0.2,
                        ID=c("SYMBOLs","ENTREZIDs","ENSEMBLIDs"),
                        colors=NULL, outdir=NULL, basename=NULL,
                        width=10, height=10
                        ){
  #Working with dataframes
  require(tidyverse)
  #STRING database
  require(STRINGdb)
  #Working with network objects
  require(igraph) #vertex_attr()
  #Graphing networks
  require(ggraph)
  require(scatterpie) #geom_scatterpie()
  require(ggnetwork) #theme_blank()
  require(scales)
  set.seed(8434)
  
  #### Data ####
  #STRING database
  message(paste("Mapping genes to STRING database version", version))
  string_db <- STRINGdb$new(version=as.character(version), species=9606,
                            score_threshold=score_threshold, input_directory="")
  
  #Format gene vector to matrix
  genes.mat <- as.matrix(genes)
  colnames(genes.mat) <- "gene"
  
  #Map genes to STRING
  map <- string_db$map(genes.mat, "gene", removeUnmappedRows = TRUE)
  
  #### Add color groups ####
  if(!is.null(enrichment)){
    message("\nFormatting enrichment colors")
    #Get significant enrichments
    col.mat <- enrichment %>% 
      dplyr::select(Description, size.overlap.term, p.adjust,
                    all_of(ID)) %>% 
      filter(size.overlap.term >= size.overlap.term & p.adjust <= p.adjust)
      
    #Error if no terms to plot
    if(nrow(col.mat)==0) {stop("No significant enrichment terms. 
                               Try increasing p.adjust.")}
    
    #Format enrichment results for scatterpie plotting
    col.mat.format <- col.mat %>% 
      #Split gene lists within terms
      mutate(gene = strsplit(as.character(get(ID)), "/")) %>% 
      unnest(gene) %>%
      #spread terms to columns
      distinct(gene, Description) %>% 
      mutate(value=1) %>% 
      #Add string ID
      left_join(map, by = "gene") %>% 
      select(-gene) %>% 
      distinct() %>% 
      #Calculate terms per ID
      group_by(STRING_id) %>% 
      mutate(total = sum(value)) %>%
      ungroup() %>% 
      #Calculate proportions within terms
      pivot_wider(names_from = Description, values_fill = 0) %>% 
      mutate(across(-c(STRING_id,total), ~./total))
    
    #Add to STRING data and create dummy group for genes without enrichment
    map.unique <- map %>% 
      #collapse gene names
      group_by(STRING_id) %>% 
      mutate(gene = paste(unique(gene), collapse=" / ")) %>% 
      #add enrichment color info
      full_join(col.mat.format, by = "STRING_id") %>%
      distinct() %>% 
      #no enrichment group
      mutate(none = ifelse(is.na(total),1,0)) %>% 
      ungroup() %>% 
      #fill in 0
      mutate(across(everything(), ~replace_na(., 0)))

  } else {
    #Collapse duplicate STRING ID
    map.unique <- map %>% 
      group_by(STRING_id) %>% 
      summarise(gene = paste(unique(gene), collapse = " / ")) %>% 
      #Dummy color group
      mutate(none=1)
  }
 
  #### Network ####
  # Create igraph object 
  subgraph <- string_db$get_subnetwork(map.unique$STRING_id)
  
  # Arrange metadata as in network
  map.arrange <- map.unique %>% 
    arrange(match(STRING_id, c(vertex_attr(subgraph)$name)))
  
  # Set attributes
  ## Check order first
  if(!identical(vertex_attr(subgraph)$name, map.arrange$STRING_id)){
    stop("igraph gene order does not match color information.")
  }
  
  ##gene names
  V(subgraph)$symbol <- map.arrange$gene
  ##enrichment colors
  for(term in colnames(map.arrange)[-c(1:3)]){
    vertex_attr(subgraph)[[term]] <- unlist(map.arrange[term])
  }

  #### Set color values ####
  if(is.null(colors)){
    color.vec <- c(hue_pal()(ncol(map.arrange)-4), "grey70")
  } else {
    color.vec <- colors
  }
  
  #### Plot ####
  message("Plotting. PLEASE IGNORE attribute warning.")
  #Get xy of nodes for manual layout
  xy <- layout_with_fr(subgraph)
  
  V(subgraph)$x <- xy[, 1]
  V(subgraph)$y <- xy[, 2]
  
  plot <- ggraph(subgraph, layout= "manual", x = V(subgraph)$x, y = V(subgraph)$y) +
    #Edges
    geom_edge_link(aes(width=combined_score), color="grey70") +
    scale_edge_width(range = c(0.2,2), name="STRING score") +
    #Nodes
    geom_scatterpie(data=as_data_frame(subgraph, "vertices"),
                    cols=sort(colnames(map.arrange)[-c(1:3)]), color=NA,
                    pie_scale = 0.7) +
    scale_fill_manual(values=color.vec,
                      name="Enrichment") +
    geom_nodetext(aes(x = V(subgraph)$x, y = V(subgraph)$y,
                      label=V(subgraph)$symbol), size=2) +
    theme_blank() + coord_fixed()
  plot(plot)
  #### Save ####
  filename <- paste(outdir,basename,"STRING.network.pdf", sep="")
  ggsave(filename, plot, height=height, width=width)
}
