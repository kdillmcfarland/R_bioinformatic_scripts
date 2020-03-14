"Boxplots of expression of a list of genes or modules

Saves individual PDF boxplots of gene/module expression by variables of
interest.

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
  voom.dat = Filepath to .csv or name of object in environment 
             containing voom normalized counts.
  pval.dat = Filepath to .csv or name of object in environment 
             containing limma results output by 'extract.pval.R'
  meta.dat = Filepath to .csv containing metadata. Only required if 
             voom.dat is csv, not voom object
  genes.toPlot = Character vector listing genes to plot
  vars = Character vector of variables in voom.dat$targets OR meta.dat
         to plot
  color.var = Variable in voom.dat$targets OR meta.dat to color points
              in plots
  outdir = Filepath to directory to save results

OPTIONAL
   interaction = Logical if should plot interaction of FIRST 2 vars.
                 Vars must be factor or character. Default is FALSE
   colors = If do not want to use default ggplot colors, list of colors
            to use. Must be same length as levels of color.var
   name = Character string to prepend to output names. Default is NULL
   gene.key = Filepath to Ensembl gene key to name genes in plots.
              Generally 'EnsemblToHGNC_GRCh38.txt'. Default is NULL
   cores = Number of parallel cores to use. Default is 1
   
Example
  plot.all(voom.dat='P259.2_voom.counts.csv', 
         pval.dat='P259.2.gene.pval.interact.csv', 
         meta.dat='P259_all_metadata.csv', 
         genes.ToPlot=c('gene1', 'gene2'),
         gene.key='EnsemblToHGNC_GRCh38.txt',
              vars=c('drug','virus'),
              interaction=TRUE,
              color.var='donorID', 
              colors=samp.cols,
              outdir='figs/gene_P259.2/', 
              name='P259.2_expression_',
              cores=3)
"

#################

plot.all <- function(voom.dat, pval.dat, meta.dat, 
                     genes.toPlot, vars,
                          interaction=FALSE,
                          color.var=NULL, colors=NULL,
                          outdir=NULL, name=NULL, 
                          gene.key=NULL,
                          cores=1){
########## SETUP ########## 
# Data manipulation and figures
library(tidyverse)
# Multi-panel figures for ggplot
library(cowplot)
#Set seed
set.seed(4389)

########## Load data ########## 
#Voom normalized counts
if(is.character(voom.dat)){
  voom.dat.loaded <- read_csv(voom.dat) %>% 
    filter(.[[1]] %in% genes.toPlot)
} else if(class(voom.dat) == "EList"){
  voom.dat.loaded <- as.data.frame(voom.dat$E) %>% 
    rownames_to_column() %>% 
    filter(rowname %in% genes.toPlot)
} else {
  stop("Voom data must be CSV on disk or EList object in environment")
}

#Pvalues
if(is.character(pval.dat)){
  pval.dat.loaded <- read_csv(pval.dat) %>% 
    dplyr::select(1, adj.P.Val, group) %>% 
    filter(.[[1]] %in% genes.toPlot)
} else if(class(pval.dat) == "data.frame"){
  pval.dat.loaded <- pval.dat %>% 
    dplyr::select(1, adj.P.Val, group)%>% 
    filter(.[[1]] %in% genes.toPlot)
} else {
  stop("P-value data must be CSV on disk or data frame in environment")
}

#Metadata
if(class(voom.dat) == "EList"){
  meta.dat.loaded <- as.data.frame(voom.dat$targets) %>% 
    dplyr::select(libID, color.var, vars)
} else if(is.character(meta.dat)){
  meta.dat.loaded <- read_csv(meta.dat) %>% 
    dplyr::select(libID, color.var, vars)
} else {
  stop("Metadata must be CSV on disk or part of EList voom object in environment.")
}

########## Format data ########## 
#Rename 1st column to match
colnames(pval.dat.loaded)[1] <- "gene"
colnames(voom.dat.loaded)[1] <- "gene"
  
# combine logCPM, meta, and pval data
plot.dat <- voom.dat.loaded %>% 
  pivot_longer(-1, names_to = "libID", 
               values_to = "voom.count") %>% 
  left_join(meta.dat.loaded, by="libID") %>% 
  left_join(pval.dat.loaded, by="gene") %>% 
  mutate(color.var = factor(get(color.var)))

########## Plots ########## 
#List all genes/modules
to_plot <- sort(unique(plot.dat$gene))

# Setup parallel computing
library(doParallel)
library(foreach)
registerDoParallel(cores=cores)

##########  Loop through genes ########## 
foreach(i = 1:length(to_plot)) %dopar% {
  print(i)
  
  #Subset data to gene/module of interest
  plot.dat.sub <- plot.dat %>% 
    filter(gene == to_plot[i])
  
########## Loop through variables ########## 
  plot_list = list()
  
  for(j in 1:length(vars)){

    #Type = factor or character variables
    if(is.factor(plot.dat.sub[[vars[j]]]) |
       is.character(plot.dat.sub[[vars[j]]])) {
      
      plot.dat.sub.fct <- plot.dat.sub %>% 
        mutate_at(vars(vars), ~fct_relevel(as.factor(.), 
                                           "none", after = 0))
      
      #List levels of variable of interest
      var.levels <- plot.dat.sub.fct %>% 
        select(vars[j]) %>% 
        distinct() %>% 
        unlist(use.names = FALSE)
      
      #Filter data to fdr values for levels assoc with variable of interest
      plot.dat.sub2 <- plot.dat.sub.fct %>% 
        filter(group %in% var.levels |
                 group == vars[j])
      
      #Extract plot title with FDR
      plot.title <- paste("FDR=", 
                          formatC(unique(plot.dat.sub2$adj.P.Val), 
                                  format = "e", digits = 4), sep="")
      
      plot1 <- plot.dat.sub2 %>% 
        ggplot(aes_string(x=vars[j], y="voom.count")) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(color=color.var), height=0, width=0.2) +
        theme_classic() +
        labs(title=plot.title, y="Normalized log2 expression",
             x="") +
        theme(legend.position = "none", 
              plot.title = element_text(size=9))
      
      if(is.null(colors)){
        plot_list[[j]] <- plot1
      } else{
        plot1 <- plot1 + scale_color_manual(values=colors)
        plot_list[[j]] <- plot1
      }
      
    } else if(is.numeric(plot.dat.sub[[vars[j]]])){
      #Filter data to fdr values for variable of interest
      plot.dat.sub2 <- plot.dat.sub %>% 
        filter(group == vars[j])
      
      #Extract plot title with FDR
      plot.title <- paste("FDR=", 
                          formatC(unique(plot.dat.sub2$adj.P.Val), 
                          format = "e", digits = 4), sep="")
      plot1 <- plot.dat.sub2 %>% 
        ggplot(aes_string(x=vars[j], y="voom.count")) +
        geom_point(aes(color=color.var)) +
        theme_classic() +
        labs(title=plot.title, y="Normalized log2 expression") +
        theme(legend.position = "none", 
              plot.title = element_text(size=9)) +
        geom_smooth(method='lm', formula= y~x, color="black",
                    se=FALSE)
      
      if(is.null(colors)){
        plot_list[[j]] <- plot1
      } else{
        plot1 <- plot1 + scale_color_manual(values=colors)
        plot_list[[j]] <- plot1
      }
    } else{
    stop("Variables of interest must be numeric, character, or factor.")
    }

  }
  
  #Interaction variable
  if(interaction){
    #List all interactions
    var.levels <- plot.dat.sub %>% 
      select(vars[1],vars[2]) %>% 
      distinct() %>% 
      mutate(interaction = paste(get(vars[1]), get(vars[2]), sep=":")) %>% 
      select(interaction) %>% 
      unlist(use.names = FALSE)
    
    #Create FDR plot title if exists in the data
    if(any(var.levels %in% plot.dat.sub$group)){
      plot.dat.sub2 <- plot.dat.sub %>% 
        filter(grepl(":", group))
      
      plot.title <- paste("FDR=", formatC(unique(plot.dat.sub2$adj.P.Val), 
                          format = "e", digits = 4), sep="")
    } else{
      plot.dat.sub2 <- plot.dat.sub %>% 
        select(-adj.P.Val,-group) %>% 
        distinct()
      
      plot.title <- " "
    }
    
    #Interaction plot
    plot2 <- plot.dat.sub2 %>% 
      ggplot(aes(x=get(vars[1]):get(vars[2]), y=voom.count)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color=color.var), height=0, width=0.2) +
      theme_classic() +
      labs(title=plot.title, y="Normalized log2 expression",
           x="") +
      theme(legend.position = "right",
            plot.title = element_text(size=9))
    
    if(is.null(colors)){
      plot2 <- plot2
    } else{
      plot2 <- plot2 + scale_color_manual(values=colors)
    }
    }
    
  #### Combine plots ####
  if (!is.null(gene.key)){
    gene.key <- read_tsv(gene.key) %>% 
      filter(ensembl_gene_id %in% plot.dat.sub$gene)
    
    plot_row0 <- ggdraw() + 
      draw_label(paste(to_plot[i], unique(gene.key$hgnc_symbol), sep=" "),
                 fontface='bold', x=0, hjust=0, vjust=4) +
      theme(plot.margin = margin(0, 0, 0, 7))
  } else{
    plot_row0 <- ggdraw() + 
      draw_label(to_plot[i],
                 fontface='bold', x=0, hjust=0, vjust=4) +
      theme(plot.margin = margin(0, 0, 0, 7))
  }
  
  plot_row1 <- plot_grid(plot_list[[1]], plot_list[[2]],
                         align="hv", nrow=1)
  plot_final <- plot_grid(plot_row0, plot_row1, plot2,
                          nrow=3, 
                          rel_heights = c(.4,1,1))
  
  #### Save to disk
  #### Include gene name if desired
  dir.create(path=outdir, showWarnings = FALSE)
  if (!is.null(gene.key)){
    filename <- paste(outdir, name,
                      unique(plot.dat.sub[,1]), "_",
                      unique(gene.key$hgnc_symbol),
                      ".pdf", sep="")
    ggsave(filename, plot_final, width=6, height=7)
  } else{
    filename <- paste(outdir, name,
                      unique(plot.dat.sub[,1]), ".pdf", sep="")
    ggsave(filename, plot_final, width=6, height=7)
  }
}

print("All plots complete.")
}