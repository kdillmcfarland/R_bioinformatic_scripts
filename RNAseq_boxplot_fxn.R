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
             voom.dat is NOT a voom object
  join.var = Variable name to use when joining data.
  genes.toPlot = Character vector listing genes to plot
  vars = Character vector of variables in voom.dat$targets OR meta.dat
         to plot
  color.var = Variable in voom.dat$targets OR meta.dat to color points
              in plots
  outdir = Filepath to directory to save results
  width, height = Dimensions of saved plot

OPTIONAL
   interaction = Logical if should plot interaction of FIRST 2 vars.
                 Vars must be factor or character. Default is FALSE
   colors = If do not want to use default ggplot colors, list of colors
            to use. Must be same length as levels of color.var
   name = Character string to prepend to output names. Default is NULL
   gene.key = Filepath to Ensembl gene key to name genes in plots.
              Generally 'EnsemblToHGNC_GRCh38.txt'. Default is NULL
   cores = Number of parallel cores to use. Default is 1
   vars_transform = Character vector of transformations to apply to x
                    variables. Can ONLY be applied to numeric variables
                    and NOT to interaction terms. Default is NULL
   
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

plot.all <- function(voom.dat, pval.dat, meta.dat=NULL, contrast.mat=NULL,
                     join.var, genes.toPlot, vars,
                          interaction=FALSE,
                          color.var=NULL, colors=NULL,
                          outdir=NULL, name=NULL, 
                          gene.key=NULL,
                          cores=1, width, height,
                          vars_transform=NULL){
########## SETUP ########## 
# Data manipulation and figures
library(tidyverse)
# Multi-panel figures for ggplot
library(cowplot)
#pval annotation
library(ggpubr)
#Set seed
set.seed(4389)

########## Load data ########## 
#Voom normalized counts
if(is.character(voom.dat)){
  voom.dat.loaded <- read_csv(voom.dat) %>% 
    filter(.[[1]] %in% genes.toPlot)
} else if(class(voom.dat) == "EList" |
          class(voom.dat) == "DGEList"){
  voom.dat.loaded <- as.data.frame(voom.dat$E) %>% 
    rownames_to_column() %>% 
    filter(rowname %in% genes.toPlot)
} else if("data.frame" %in% class(voom.dat)){
  voom.dat.loaded <- voom.mods %>% 
    rownames_to_column() %>% 
    filter(rowname %in% genes.toPlot)
} else {
  stop("Voom data must be CSV on disk or edgeR/data.frame object in environment")
}

#Pvalues
if(is.character(pval.dat)){
  pval.dat.loaded <- read_csv(pval.dat) %>% 
    dplyr::select(1, adj.P.Val, group) %>% 
    filter(.[[1]] %in% genes.toPlot)
} else if("data.frame" %in% class(pval.dat)){
  pval.dat.loaded <- pval.dat %>% 
    dplyr::select(1, adj.P.Val, group)%>% 
    filter(.[[1]] %in% genes.toPlot)
} else {
  stop("P-value data must be CSV on disk or data.frame in environment")
}

#Metadata
if(class(voom.dat) == "EList" |
   class(voom.dat) == "DGEList"){
  meta.dat.loaded <- as.data.frame(voom.dat$targets) %>% 
    dplyr::select(join.var, all_of(color.var), all_of(vars))
} else if(is.character(meta.dat)){
  meta.dat.loaded <- read_csv(meta.dat) %>% 
    dplyr::select(join.var, all_of(color.var), all_of(vars))
} else if("data.frame" %in% class(meta.dat)){
  meta.dat.loaded <- meta.dat %>% 
    dplyr::select(join.var, all_of(color.var), all_of(vars))
} else {
  stop("Metadata must be CSV on disk or part of EList/data.frame object in environment.")
}

########## Format data ########## 
#Rename 1st column to match
colnames(pval.dat.loaded)[1] <- "gene"
colnames(voom.dat.loaded)[1] <- "gene"
  
# combine logCPM, meta, and pval data
plot.dat <- voom.dat.loaded %>% 
  pivot_longer(-1, names_to = join.var, 
               values_to = "voom.count") %>% 
  left_join(meta.dat.loaded) %>% 
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
foreach(i = 1:length(to_plot), .verbose = TRUE) %dopar% {
  
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
        mutate_at(vars(vars[j]), ~fct_relevel(as.factor(.), 
                                           "none", after = 0))
      
      #List levels of variable of interest
      var.levels <- plot.dat.sub.fct %>% 
        select(vars[j]) %>% 
        distinct() %>% 
        unlist(use.names = FALSE)
      
      #IF variable was in model
      if(vars[j] %in% plot.dat.sub.fct$group){
      #Filter data to fdr values for levels assoc with variable of interest
      plot.dat.sub2 <- plot.dat.sub.fct %>% 
        filter(group %in% var.levels |
                 group == vars[j])
      #Extract plot title with FDR
      plot.title <- paste("FDR=", 
                          formatC(unique(plot.dat.sub2$adj.P.Val), 
                                  format = "e", digits = 4), sep="")
      } else{
        plot.dat.sub2 <- plot.dat.sub.fct %>% 
          select(-group, -adj.P.Val) %>% 
          distinct()
        #Make plot title w/o FDR
        plot.title <- "Not in model"
      }
      
      plot1 <- plot.dat.sub2 %>% 
        ggplot(aes_string(x=vars[j], y="voom.count")) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(color=color.var), height=0, width=0.2) +
        theme_classic() +
        labs(title=plot.title, y="Normalized log2 expression") +
        theme(legend.position = "none", 
              plot.title = element_text(size=9))
      
      if(is.null(colors)){
        plot_list[[j]] <- plot1
      } else{
        plot1 <- plot1 + scale_color_manual(values=colors)
        plot_list[[j]] <- plot1
      }
      
    } else 
    #Type = numeric
    if(is.numeric(plot.dat.sub[[vars[j]]])){
      #IF variable was in model
      if(vars[j] %in% plot.dat.sub$group){
        #Filter data to fdr values for variable
        plot.dat.sub2 <- plot.dat.sub %>% 
          filter(group == vars[j])
        #Extract plot title with FDR
        plot.title <- paste("FDR=", 
                            formatC(unique(plot.dat.sub2$adj.P.Val), 
                                    format = "e", digits = 4), sep="")
        
        ##### Transform if provided
        if(!is.null(vars_transform)){
          transformation <- gsub("x", vars[j], vars_transform[j])
          
          plot.dat.sub2 <- plot.dat.sub2 %>% 
            #transform
            mutate_(transformation = transformation)
          x.var <- transformation
          
        } else {
          x.var <- vars[j]
        }
      } 
      #Else variable was not in model
      else{
        plot.dat.sub2 <- plot.dat.sub %>% 
          select(-group, -adj.P.Val) %>% 
          distinct()
        #Make plot title w/o FDR
        plot.title <- "Not in model"
        
        ##### Transform if provided
        if(!is.null(vars_transform)){
          transformation <- gsub("x", vars[j], vars_transform[j])
          
          plot.dat.sub2 <- plot.dat.sub2 %>% 
            #transform
            mutate_(transformation = transformation)
          x.var <- transformation
          
        } else {
          x.var <- vars[j]
        }
      }
       
        plot1 <- plot.dat.sub2 %>% 
        ggplot(aes_string(x=x.var, y="voom.count")) +
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
    
    var.levels.addtl <- c(paste(vars[1],vars[2],sep=":"),
                          paste(vars[2],vars[1],sep=":"))
    
    if(!is.null(contrast.mat)){
      contrast.levels <- colnames(contrast.mat)
      var.levels.all <- c(var.levels, var.levels.addtl, contrast.levels)
    } else {
      var.levels.all <- c(var.levels, var.levels.addtl)
    }
    
    #Create FDR plot title if exists in the data
    if(any(var.levels.all %in% plot.dat.sub$group)){
      title.dat <- plot.dat.sub %>% 
        filter(group %in% var.levels.all) %>% 
        select(group, adj.P.Val) %>% 
        distinct()

      plot.title <- paste(title.dat$group, "FDR =", 
                          formatC(title.dat$adj.P.Val,
                                  format = "e", digits = 4), 
                          sep=" ",collapse="\n")
      
      plot.dat.sub2 <- plot.dat.sub %>% 
        filter(group %in% var.levels.all) %>% 
        select(-adj.P.Val,-group) %>% 
        distinct()
      
    } else{
      plot.dat.sub2 <- plot.dat.sub %>% 
        select(-adj.P.Val,-group) %>% 
        distinct()
      
      plot.title <- " "
    }
    
    #Interaction plot
    plot2 <- plot.dat.sub2 %>% 
      ggplot(aes(x=paste(get(vars[1]),get(vars[2]), sep=":"), 
                 y=voom.count)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color=color.var), height=0, width=0.2) +
      theme_classic() +
      labs(title=plot.title, y="Normalized log2 expression",
           x="", color=color.var) +
      theme(legend.position = "right",
            plot.title = element_text(size=9))

##### Color points #####   
    if(!is.null(colors)){
      plot2 <- plot2 + 
        scale_color_manual(values=colors)
    }
    
#### Add signif bars for contrasts #####     
    if(!is.null(contrast.mat)){
      
      stat.dat <- title.dat %>% 
        filter(adj.P.Val <= 0.05) 
      
      if(nrow(stat.dat)>0){
        
      stat.dat <- stat.dat %>% 
        separate(group, into=c("group1","group2"), sep=" - ") %>% 
        mutate(group1=gsub("_|-", ":", group1),
               group2=gsub("_|-", ":", group2)) %>% 
        mutate(adj.P.Val = formatC(adj.P.Val, 
                                   format = "e", digits = 4))
      
      y.pos <- c() 
      for(k in 1:nrow(stat.dat)){
        scale.val <- k/10+1
        temp <- max(plot.dat.sub2[,"voom.count"])*scale.val
        
        y.pos <- c(y.pos, temp)
      }
      
      stat.dat <- stat.dat %>% 
        mutate(y.position = y.pos)
        
      plot2 <- plot2 +
        stat_pvalue_manual(stat.dat, label="adj.P.Val")
    }}
    
  } else{
      plot2 <- NULL
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
  
  plot_row1 <- plot_grid(plotlist = plot_list,
                         align="hv", nrow=1)
  
  if(!is.null(plot2)){
    plot_final <- plot_grid(plot_row0, plot_row1, plot2,
                            nrow=3, 
                            rel_heights = c(.4,1,1))
  } else{
    plot_final <- plot_grid(plot_row0, plot_row1,
                            nrow=2, 
                            rel_heights = c(.4,1))
  }
  
  #### Save to disk
  #### Include gene name if desired
  dir.create(path=outdir, showWarnings = FALSE)
  if (!is.null(gene.key)){
    filename <- paste(outdir, name,
                      unique(plot.dat.sub[,1]), "_",
                      unique(gene.key$hgnc_symbol),
                      ".pdf", sep="")
    ggsave(filename, plot_final, width=width, height=height)
  } else{
    filename <- paste(outdir, name,
                      unique(plot.dat.sub[,1]), ".pdf", sep="")
    ggsave(filename, plot_final, width=width, height=height)
  }
}

print("All plots complete.")
}