"Create STRING network of gene list

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
  layout = layout type of network. One of Fruchterman-Reingold 'fr', bipartite 'bipar',
           star, tree, circle, Kamada-Kawai 'kk', graphopt, gem, Davidson-Harel 'dh',
           sphere, grid, large 'lgl', mds, Sugiyama 'sugi' Default is Fruchterman-Reingold
OPTIONAL
Coloring by enrichment terms
  enrichment = Matrix/data frame output by hypergeo_enricher function. Necessary
               columns are Description, size.overlap.term, p.adjust, and IDs (see below)
  discard = Character vector listing groups of nodes to remove from network. Can be 
            'edge.keep.enrich' to remove nodes without any edge connections but keep
            unconnected nodes with enrichment or 'edge' to remove nodes without edges
            regardless of enrichment. Default is 'none'
  ID = Character for ID type of gene list. One of c('SYMBOLs','ENTREZIDs','ENSEMBLIDs')
  size.overlap.term = Gene overlap minimum for enrichment terms to include in plot.
                      Default is 2 
  p.adjust = FDR maximum for enrichment terms to include in plot. Default is 0.2
  colors = Vector of colors for manual scales. Default it to use ggplot colors.
  node.text.size = Numeric
Saving plot
  basename = Character to prepend to output names.
  outdir = Character for output directory for plot
  height, width = Plot dimensions in inches
"

#################
string.plot <- function(genes, version="11", score_threshold=400,
                        layout='fr',
                        enrichment=NULL, size.overlap.term=2, p.adjust=0.2,
                        discard="none",
                        ID=c("SYMBOLs","ENTREZIDs","ENSEMBLIDs"),
                        colors=NULL, node.text.size=2, pie_scale=1,
                        outdir=NULL, basename=NULL,
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
      ungroup() %>% 
      dplyr::filter(size.overlap.term >= size.overlap.term & p.adjust <= p.adjust) %>% 
      dplyr::select(Description, all_of(ID))
      
    #Error if no terms to plot
    if(nrow(col.mat)==0) {stop("No significant enrichment terms. 
                               Try increasing p.adjust.")}
    
    #Format enrichment results for scatterpie plotting
    col.mat.format <- col.mat %>% 
      #Split gene lists within terms
      dplyr::mutate(gene = strsplit(as.character(get(ID)), "/")) %>% 
      tidyr::unnest(gene) %>%
      #spread terms to columns
      distinct(gene, Description) %>% 
      dplyr::mutate(value=1) %>% 
      #Add string ID
      inner_join(map, by = "gene") %>% 
      dplyr::select(-gene) %>% 
      dplyr::distinct() %>% 
      #Calculate terms per ID
      group_by(STRING_id) %>% 
      dplyr::mutate(total = sum(value)) %>%
      ungroup() %>% 
      #Calculate proportions within terms
      arrange(match(Description,unique(enrichment$Description))) %>% 
      pivot_wider(names_from = Description, values_fill = 0) %>% 
      dplyr::mutate(across(-c(STRING_id,total), ~./total))
    
    #Add to STRING data and create dummy group for genes without enrichment
    map.unique <- map %>% 
      #collapse gene names
      group_by(STRING_id) %>% 
      dplyr::mutate(gene = paste(unique(gene), collapse=" / ")) %>% 
      #add enrichment color info
      full_join(col.mat.format, by = "STRING_id") %>%
      distinct() %>% 
      #no enrichment group
      dplyr::mutate(none = ifelse(is.na(total),1,0)) %>% 
      ungroup() %>% 
      #fill in 0
      dplyr::mutate_if(is.numeric, ~replace_na(., 0))

  } else {
    #Collapse duplicate STRING ID
    map.unique <- map %>% 
      group_by(STRING_id) %>% 
      dplyr::mutate_if(is.numeric, ~replace_na(., 0))
  }
 
  #### Network ####
  # Create igraph object 
  subgraph <- string_db$get_subnetwork(map.unique$STRING_id)
  
  #Discard unconnected and uncolored if specified
  ##All nodes
  nodes <- which(degree(subgraph)>=0)
  ##Nodes without edge connections
  isolated <- which(degree(subgraph)==0)
  

  ##Remove unconnected regardless of enrichment
  if(discard == "edge"){
    subgraph.filter <- delete.vertices(subgraph, isolated)
  } else if(discard == "edge.keep.enrich"){
    ##Nodes without enrichment
    unenrich <- map.unique %>% 
      dplyr::filter(none==1) %>% 
      distinct(STRING_id) %>% unlist(use.names = FALSE)
  ##Remove unconnected that are also unenriched
    isolated.unenrich <- isolated[names(isolated) %in% unenrich]
    subgraph.filter <- delete.vertices(subgraph, isolated.unenrich)
  } else if(discard == "none") {
    subgraph.filter <- subgraph
  }
  
  # Arrange metadata as in network
  map.arrange <- map.unique %>% 
    dplyr::filter(STRING_id %in% vertex_attr(subgraph.filter)$name) %>% 
    arrange(match(STRING_id, c(vertex_attr(subgraph.filter)$name)))

  # Set attributes
  ## Check order first
  if(!identical(vertex_attr(subgraph.filter)$name, map.arrange$STRING_id)){
    stop("igraph gene order does not match color information.")
  }
  
  ##gene names
  V(subgraph.filter)$symbol <- map.arrange$gene
  ##enrichment colors
  for(term in colnames(map.arrange)[-c(1:3)]){
    vertex_attr(subgraph.filter)[[term]] <- unlist(map.arrange[term])
  }

  #### Set color values ####
  if(!is.null(colors)){
    color.vec <- colors
  } else if(!is.null(enrichment)){
    if("none" %in% colnames(col.mat.format)){
      #ggplot colors for sets in plot
      color <- hue_pal()(ncol(col.mat.format)-2)
      #Replace none with grey
      all.term <- sort(colnames(col.mat.format)[-c(1:2)])
      none.index <- match("none", all.term)
      color.vec <- c(color[1:none.index-1], "grey70",
                     color[none.index:length(color)])
    } else {
      color.vec <- hue_pal()(ncol(col.mat.format)-2)
    }
  } else if (is.null(enrichment) & is.null(colors)){
    color.vec <- "#d9d9d9"
  }
  
  #### Plot ####
  message("\nPlotting. PLEASE IGNORE attribute warning.")
  #Get xy of nodes for manual layout
  set.seed(8434)
  ##set layout
  if(layout == "fr"){ xy <- layout_with_fr(subgraph.filter) } else
    if(layout == "bipar"){ xy <- layout_as_bipartite(subgraph.filter) } else
      if(layout == "star"){ xy <- layout_as_star(subgraph.filter) } else
        if(layout == "tree"){ xy <- layout_as_tree(subgraph.filter) } else
          if(layout == "circle"){ xy <- layout_in_circle(subgraph.filter) } else
            if(layout == "kk"){ xy <- layout_with_kk(subgraph.filter) } else
              if(layout == "graphopt"){ xy <- layout_with_graphopt(subgraph.filter) } else
                if(layout == "gem"){ xy <- layout_with_gem(subgraph.filter) } else
                  if(layout == "dh"){ xy <- layout_with_dh(subgraph.filter) } else
                    if(layout == "sphere"){ xy <- layout_on_sphere(subgraph.filter) } else
                      if(layout == "grid"){ xy <- layout_on_grid(subgraph.filter) } else
                        if(layout == "lgl"){ xy <- layout_with_lgl(subgraph.filter) } else
                          if(layout == "mds"){ xy <- layout_with_mds(subgraph.filter) } else
                            if(layout == "sugi"){ xy <- layout_with_sugiyama(subgraph.filter) }

  V(subgraph.filter)$x <- xy[, 1]
  V(subgraph.filter)$y <- xy[, 2]
  
  plot <- ggraph(subgraph.filter, layout= "manual", 
                 x = V(subgraph.filter)$x, y = V(subgraph.filter)$y) +
    #Edges
    geom_edge_link(aes(width=combined_score), color="grey70") +
    scale_edge_width(range = c(0.2,2), name="STRING score") 
  
  #Add nodes
  if(!is.null(enrichment)){
    plot.col <- plot + 
      geom_scatterpie(data=as_data_frame(subgraph.filter, "vertices"),
                      cols=colnames(map.arrange)[-c(1:3)], color=NA,
                      pie_scale = pie_scale) +
      scale_fill_manual(values=color.vec, name="Enrichment") +
      geom_nodetext(aes(x = V(subgraph.filter)$x, y = V(subgraph.filter)$y,
                        label=V(subgraph.filter)$symbol), 
                    size=node.text.size) +
      theme_blank() + coord_fixed()
  } else{
    plot.col <- plot + 
      geom_nodes(aes(x = V(subgraph.filter)$x, y = V(subgraph.filter)$y,
                     fill=NULL), size = pie_scale*10, color=color.vec) +
      geom_nodetext(aes(x = V(subgraph.filter)$x, y = V(subgraph.filter)$y,
                          label=V(subgraph.filter)$symbol), 
                    size=node.text.size) +
      theme_blank() + coord_fixed()
  }
  
  plot(plot.col)
  
  #### Save ####
  filename <- paste(outdir,basename,"STRING.network.pdf", sep="")
  ggsave(filename, plot.col, height=height, width=width)
}
