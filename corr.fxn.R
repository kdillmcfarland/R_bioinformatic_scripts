################
"
Run correlation of counts data and associated metadata

Automatically subsets to complete samples.

REQUIRED
  counts = counts matrix or data frame with genes/modules as columns and
           samples as rows
  metadata = metadata for samples in count data
  match.by = column name from metadata that matches sample names in counts.
             Default it 'libID'
  corr = Correlation metric for Hmisc::rcorr( ). Default is 'pearson'
  basename = Prefix for results naming. Default is 'counts_meta'

OPTIONAL (Default is NULL / FALSE)
  subset_count = If need to calculate correlations for subset of 
                 genes/modules, provide their row names
  subset_meta = If need to calculate correlations for subset of metadata,
                 provide their column names
  plot = Logical if create heatmaps. Default is FALSE
  heatmap_meta = Additional metadata to use to annotation heatmaps. Must have
              first column values to match count column names
                 
EXAMPLE
corr.fxn(counts = counts, metadata = meta, match.by='libID',
         corr='Pearson')
"
################

corr.fxn <- function(counts, metadata, match.by="libID", corr="pearson", 
                     subset_count=NULL, subset_meta=NULL, 
                     basename="counts_meta", 
                     heatmap=FALSE, heatmap_meta=NULL,
                     dotplot=FALSE, color.visit=FALSE){
  require(tidyverse, quietly = TRUE)
  require(Hmisc, quietly = TRUE)
  require(corrplot, quietly = TRUE)
  require(plotly, quietly = TRUE)
  require(heatmaply, quietly = TRUE)
  require(ComplexHeatmap, quietly = TRUE)
  require(circlize, quietly = TRUE)
  set.seed(4389)

##### Filter to genes/modules of interest #####
  if(!is.null(subset_count)){
    count.dat <- counts %>% 
      dplyr::select(1, all_of(subset_count))
    
  } else {
    count.dat <- counts
  }

##### Filter to plot metadata if available #####
  if(!is.null(heatmap_meta)){
    heatmap_meta.dat <- heatmap_meta %>% 
      rename("rowname" = 1) %>% 
      arrange(rowname) %>% 
      filter(rowname %in% colnames(count.dat))
  }
  
##### Filter to metadata of interest #####
  if(!is.null(subset_meta)){
    meta.dat <- metadata %>% 
      dplyr::select(1, all_of(subset_meta))
  } else {
    meta.dat <- metadata
  } 
##### Filter to samples with both count and meta data#####
  shared.samp <- dplyr::intersect(count.dat[,match.by], 
                                  meta.dat[,match.by]) %>% 
    unlist(use.names = FALSE)
  
  count.sub <- count.dat %>% 
    filter_at(match.by, any_vars(. %in% shared.samp)) %>% 
    arrange(get(match.by)) %>% 
    column_to_rownames(match.by) %>% 
    as.matrix()
  meta.sub <- meta.dat %>% 
    filter_at(match.by, any_vars(. %in% shared.samp)) %>% 
    arrange(get(match.by)) %>% 
    column_to_rownames(match.by) %>% 
    as.matrix()
  
print(paste("Correlating", length(shared.samp), "samples", sep=" "))

##### Correlation #####

corr.result <- rcorr(count.sub, meta.sub, type=corr)

corr.R <- as.data.frame(corr.result$r) %>% 
  rownames_to_column() %>% 
  filter(rowname %in% colnames(count.sub)) %>% 
  dplyr::select(rowname, all_of(colnames(meta.sub))) %>% 
  column_to_rownames("rowname") %>% 
  as.matrix()

corr.P <- as.data.frame(corr.result$P) %>% 
  rownames_to_column() %>% 
  filter(rowname %in% colnames(count.sub)) %>% 
  dplyr::select(rowname, all_of(colnames(meta.sub))) %>% 
  column_to_rownames("rowname") %>% 
  as.matrix()

##### Save #####
  dir.create("results/correlation", showWarnings = FALSE)
  base.filename <- paste(corr, "_", basename, ".csv", sep="")

  as.data.frame(corr.R) %>% 
    rownames_to_column("variable") %>% 
    write_csv(path=paste("results/correlation/", "R_", base.filename, sep=""))

  as.data.frame(corr.P) %>% 
    rownames_to_column("variable") %>% 
    write_csv(path=paste("results/correlation/", "P_", base.filename, sep=""))
  
##### Plots #####
  if(heatmap){
    fig.dir <- paste("figs/correlation/", basename, sep="")
    dir.create(fig.dir, showWarnings = FALSE)
    
    #### HTML heatmap ####
    print("HTML heatmap")
    html.name <- paste(fig.dir, "/all_", corr, "_",
                       basename, ".html", sep="")
    # Add Pvalue hover test
    hover.text <- as.data.frame(corr.P) %>% 
      rownames_to_column() %>% 
      mutate_at(vars(-rowname), 
                ~paste("P = ", round(., digits=4), sep=""))
    
    #Add additional annotations to hover text if provided
    if(!is.null(heatmap_meta)){
      #Select categorical variables for annotation
      annot_ly <- heatmap_meta.dat %>% 
        select_if(negate(is.numeric)) %>% 
        column_to_rownames()
      
      #plot
      heatmaply_cor(corr.R, custom_hovertext = hover.text,
                    row_side_colors = annot_ly,
                    plot_method = "ggplot",
                    file=html.name,
                    hclust_method="average")
      
    } else{
      heatmaply_cor(corr.R, custom_hovertext = hover.text,
                    plot_method = "ggplot",
                    file=html.name,
                    hclust_method="average")
    }
    #Remove heatmap temp files
    unlink(list.files(path=fig.dir, pattern="_files", 
                      full.names = TRUE), 
           recursive = TRUE)
    
    #### Static heatmap 1 ####
    print("Static heatmaps")
    pdf.name1 <- paste(fig.dir, "/all_", corr, "_",
                      basename, "1.pdf", sep="")
    set.width <- ifelse(ncol(corr.P)>6, ncol(corr.P), 6)
    
    heatmap_meta.dat.rows <- heatmap_meta.dat %>% 
      dplyr::select(-adj.P.Val) %>% 
      column_to_rownames()
    
    annot_rows <- rowAnnotation(df=heatmap_meta.dat.rows,
                col = list("Fold change" =
                c("Significant decrease in expression with visit"="darkblue",
                  "Not significant for visit"="darkgrey",
                  "Significant increase in expression with visit"="darkred"),
                "cell type" = c("EOS"="#F8766D",
                                "LYM"="#A3A500",
                                "MONO"="#00BF7D",
                                "PMN"="#00B0F6",
                                "unassign"="#E76BF3"),
                "genes in module" = colorRamp2(c(0, 1000), 
                                                c("white", "black"))))
    
    pdf(file = pdf.name1, height=nrow(corr.P)/2, width=set.width)
    draw(Heatmap(as.matrix(corr.R),
            name = corr,
            right_annotation = annot_rows))
    dev.off()
    
    #### Static heatmap 2 ####
    pdf.name2 <- paste(fig.dir, "/all_", corr, "_",
                       basename, "2.pdf", sep="")
    
    pdf(file = pdf.name2, height=nrow(corr.P)/2, width=set.width)
    
    corrplot(as.matrix(corr.R), method="color", order="original",
             #Reverse color order from default
             col=colorRampPalette(c("darkblue", "white",
                                    "darkred"))(20),
             #Change labels
             tl.col="black", tl.srt=90, cl.lim=c(-1,1),
             #Change colorlegend
             cl.pos="b", cl.length = 5,
             #Add significant labels
             p.mat = as.matrix(corr.P), insig = "label_sig",
             sig.level = c(0.01, 0.05, 0.1), 
             pch.cex=0.9, pch.col="white",
             #Add title
             title="P *<0.1 **<0.05 ***<0.01", 
             mar = c(0,0,2,0))
    
    dev.off()
  }
  if(dotplot){
    fig.dir <- paste("figs/correlation/", basename, sep="")
    dir.create(fig.dir, showWarnings = FALSE)   
    #### dot plots ####
    print("Dot plots")
    plot_dat <- as.data.frame(count.sub) %>% 
      rownames_to_column() %>% 
      pivot_longer(-rowname, names_to = "name",values_to = "count.value")
    
    plot_dat <- as.data.frame(meta.sub) %>% 
      rownames_to_column() %>%
      pivot_longer(-rowname, names_to = "meta", values_to = "meta.value") %>% 
      inner_join(plot_dat, by="rowname")
    
    for (meta.var in unique(colnames(corr.R))){
      print(meta.var)
      for (count.var in unique(rownames(corr.R))){
        
        #Subset data to modules of interest
        plot_dat.sub <- plot_dat %>% 
          filter(meta == meta.var & name == count.var)
        
        #Group labels with corr value and P
        corr.R.sub <- as.data.frame(corr.R) %>% 
          rownames_to_column() %>% 
          filter(rowname == count.var) %>% 
          dplyr::select(rowname, all_of(meta.var))
        
        corr.P.sub <- as.data.frame(corr.P) %>% 
          rownames_to_column() %>% 
          filter(rowname == count.var) %>% 
          dplyr::select(rowname, all_of(meta.var))
        
        plot.title <- paste(corr, " R = ", round(corr.R.sub[,2],
                                                  digits=3),
                            "\nP = ", round(corr.P.sub[,2], digits=3),
                            sep="")
        
        if(color.visit){
          plot1 <- plot_dat.sub %>% 
            separate(rowname, into=c("donorID", "visit"), sep="_") %>% 
            
            ggplot(aes(x=meta.value, y=count.value)) +
            #Indiv patient points
            geom_point(aes(color=visit)) +
            #Overall black trend line
            geom_smooth(formula=y~x,
                        color="black", method="lm", se=FALSE) +
            #Beautify
            theme_classic() +
            labs(x=unique(plot_dat.sub$meta), 
                 y=unique(plot_dat.sub$name), 
                 title=plot.title)
        } else {
          plot1 <- ggplot(plot_dat.sub, 
                          aes(x=meta.value, y=count.value)) +
            #Indiv patient points
            geom_point(aes(color=rowname)) +
            #Overall black trend line
            geom_smooth(formula=y~x,
                        color="black", method="lm", se=FALSE) +
            #Beautify
            theme_classic() +
            labs(x=unique(plot_dat.sub$meta), 
                 y=unique(plot_dat.sub$name), 
                 title=plot.title,
                 color=match.by)
        } 
        
        plot.name <- paste(fig.dir, "/",
                           unique(plot_dat.sub$meta), "_",
                           unique(plot_dat.sub$name), 
                           ".pdf", sep="")
        ggsave(plot.name, plot1, height=5, width=5)
      }}}
}

##### FIN #####