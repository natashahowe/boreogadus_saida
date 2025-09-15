# Neutral & Adaptive PCA & Admixture Plots Combined

# ADMIXTURE ANALYSIS

library(tidyverse)
library(reshape2)
library(patchwork)
library(tibble)
library(here) 
library(viridis)
library(stringr)
library(tools)

## Filepaths ################################################################
BAMFILE <- "./data/R/boreogadus_filtered_downsampled_bams.txt"
METADATAFILE <- "./data/R/boreogadus_metadata.csv"
COLORFILE <- paste0(here(),"./data/R/color_metadata_allpops_flip.txt")

### Metadata ###############################################################
meta_df <- read.csv(METADATAFILE, header = T) %>%
  mutate(Population = ifelse(Population=='Southhampton Island', 'Southampton Island', Population),
         Pop = case_when(Population == 'SE Baffin' ~ 'SEBaffin',
                         TRUE ~ gsub(' .*$','',Population))) %>%
  filter(Pop != 'Unknown')

color_df <- read.delim2(COLORFILE, header = T, row.names = NULL, sep = "\t") %>%
  mutate(Pop = gsub('hh','h',Pop)) %>%
  left_join(meta_df %>% distinct(Population, Pop), by = "Pop")

# call in bams
bam_df <- read.table(BAMFILE, header = F) %>%
  mutate(sampleID = basename(V1) %>% tools::file_path_sans_ext(.),
         sampleID = gsub("_sorted",'',sampleID) %>% gsub("_downsampled",'',.),
         sampleID = gsub("boreogadus_",'',sampleID)) %>%
  select(sampleID)

# call in some metadata
pop_df <- meta_df %>% # remove everything after space
  left_join(bam_df, ., by = "sampleID") 

# create mypalette
mypalette <- as.vector(color_df$Color2) # turn colors into vector
  names(mypalette) <- color_df$Pop # attach pop name to palette color
  pop.factor.levels <- color_df$Pop

## ADMIX ##############################################################
 
# before
# admix_colors <- c("K1" = "#000004","K2" = "#FCFDBF","K3" = "#DE4968",
#                   "K4" = "#8C2981","K5" = "#FE9F6D","K6" = "#3B0F70")

admix_colors <- c('navy','skyblue2','#2C967D','violetred4','#FDAA5C','khaki2')

### Neutral Admix #####################################################
admix_neutral_plot <- list()

for(K in 2:4){
  
  WILDCARD = paste0('*unfused_neutral_alpha0.05*k',K,"-1*.qopt")
  ADMIXFILE <- Sys.glob(file.path(here(),"results","admix",WILDCARD))
  basename(ADMIXFILE)
  
  admix_ks <- read.delim2(ADMIXFILE, sep = "", row.names = NULL, header = F,
                          col.names = paste0('K',1:K))
  
  admix_df <- as.data.frame(lapply(admix_ks, as.numeric))
  
  #### Prepare data for plotting ############
  
  # combine FID and admix output file to denote indiv and pop
  admix_df = as_tibble(cbind(pop_df, admix_df))
  admix_df$Pop <- factor(admix_df$Pop, levels = unique(color_df$Pop))
  
  # make column names say K instead of V
  admix_col_names <- c("Indiv", "Pop_Old", "Population", "Pop", "K1")
  # account for number of Ks in column naming
  for(i in 2:K){
    admix_col_names <- append(admix_col_names, paste0("K",i))
  }
  colnames(admix_df) <- admix_col_names
  
  # change structure of admix_df to plot
  admixed <- admix_df %>%
    select(-Pop_Old) %>%
    melt(., id = c("Indiv", "Population", "Pop"))
  colnames(admixed) <- c("Indiv", "Population", "Pop", "Cluster", "Proportion")
  
  # arrange individuals in population by their max cluster for better viewing
  order_by_cluster <- admixed %>%
    group_by(Indiv, Pop) %>%
    summarise(
      dom_Cluster = Cluster[which.max(Proportion)],
      max_prop = max(Proportion),
      .groups = "drop"
    ) %>%
    group_by(Pop) %>%
    arrange(Pop, dom_Cluster, desc(max_prop)) %>%
    mutate(within_pop_order = row_number()) %>%
    ungroup()
  
  admixed <- admixed %>%
    left_join(order_by_cluster %>% select(Indiv, within_pop_order), by = "Indiv")
  
  pop_levels <- unique(admixed$Pop[order(admixed$Pop)])
  
  admixed <- admixed %>%
    arrange(factor(Pop, levels = pop_levels), within_pop_order) %>%
    mutate(Indiv = factor(Indiv, levels = unique(Indiv)))
  
  #### Plot #################################################
  # population bar

  ## add a column called Locality_text that just has Locality as all the values (x-axis)
  admixed <- admixed %>% 
    mutate(Locality_text = "Locality")
  
  # admix data
  admix_neutral_plot[[K]] <- ggplot(admixed, aes(x = Indiv, y = Proportion, fill = Cluster)) +
    geom_bar(width = 1, stat = "identity", position = "stack") +
    facet_grid( ~ Pop, scales = "free_x", space = "free_x") +
    labs(y = paste0("K",K)) +
    theme_minimal() +
    theme(axis.line = element_blank(), 
          title = element_blank(),
          legend.text = element_text(size = 16), 
          legend.title = element_text(size = 18),
          plot.margin = unit(c(t = 0,r = 0,l = 0,b = 0), "cm"),
          axis.text.y = element_text(size = 12, color = 'black'),
          axis.text.x = element_blank(), axis.ticks.y = element_line(),
          axis.title.y = element_text(hjust=1, angle=0, size = 18),
          axis.title.x = element_blank(), strip.text.x = element_blank()) +
    scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1),
                       expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = admix_colors)
  
}

# plot a bar that represents population ids for the individuals in rows of the genotype heatmap
pop_plot <- ggplot(admixed, aes(x = Indiv, y = Locality_text, fill = Pop)) + 
  geom_tile() + 
  facet_grid( ~ Pop, scales = "free_x", space = "free_x", switch = "x") +
  labs(x = "Population", fill = "Population") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        title = element_blank(),
        strip.placement = "inside", strip.background.x = element_blank(),
        legend.position = "right", 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 18),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = mypalette)  

admix_neutral_plot[-length(admix_neutral_plot)] <- lapply(admix_neutral_plot[-length(admix_neutral_plot)], 
                                          function(p) p + theme(legend.position = "none"))

neutral_admix_K23 <- ((admix_neutral_plot[[2]] / admix_neutral_plot[[3]]) / pop_plot) + 
  plot_layout(guides = "collect",heights = c(1,1,0.25)) +
  plot_annotation(title = 'Unfused: Neutral Loci',
                  theme = theme(plot.title = element_text(hjust=0.1, size=20,
                                                          margin = margin(0,0,10,0))))
neutral_admix_K23 

neutral_admix_K24 <- Reduce('/',Filter(Negate(is.null), admix_neutral_plot)) / pop_plot + 
  plot_layout(guides = "collect",
              heights = c(rep(1,length(Filter(Negate(is.null), admix_neutral_plot))),0.25)) +
  plot_annotation(title = 'Unfused: Neutral Loci',
                  theme = theme(plot.title = element_text(hjust=0.1, size=20,
                                                          margin = margin(0,0,10,0))))
neutral_admix_K24 

## Adaptive Admix ##############################################################

admix_adaptive_plot <- list()

for(K in 2:4){
  
  WILDCARD = paste0('*_fused_adaptive_alpha0.05*k',K,"-1*.qopt")
  ADMIXFILE <- Sys.glob(file.path(here(),"results","admix",WILDCARD))
  basename(ADMIXFILE)
  
  admix_ks <- read.delim2(ADMIXFILE, sep = "", row.names = NULL, header = F,
                          col.names = paste0('K',1:K))
  
  admix_df <- as.data.frame(lapply(admix_ks, as.numeric))
  
  #### Prepare data for plotting ############
  
  # combine FID and admix output file to denote indiv and pop
  admix_df = as_tibble(cbind(pop_df, admix_df))
  admix_df$Pop <- factor(admix_df$Pop, levels = unique(color_df$Pop))
  
  # make column names say K instead of V
  admix_col_names <- c("Indiv", "Pop_Old", "Population", "Pop", "K1")
  # account for number of Ks in column naming
  for(i in 2:K){
    admix_col_names <- append(admix_col_names, paste0("K",i))
  }
  colnames(admix_df) <- admix_col_names
  
  # change structure of admix_df to plot
  admixed <- admix_df %>%
    select(-Pop_Old) %>%
    melt(., id = c("Indiv", "Population", "Pop"))
  colnames(admixed) <- c("Indiv", "Population", "Pop", "Cluster", "Proportion")
  
  # arrange individuals in population by their max cluster for better viewing
  order_by_cluster <- admixed %>%
    group_by(Indiv, Pop) %>%
    summarise(
      dom_Cluster = Cluster[which.max(Proportion)],
      max_prop = max(Proportion),
      .groups = "drop"
    ) %>%
    group_by(Pop) %>%
    arrange(Pop, dom_Cluster, desc(max_prop)) %>%
    mutate(within_pop_order = row_number()) %>%
    ungroup()
  
  admixed <- admixed %>%
    left_join(order_by_cluster %>% select(Indiv, within_pop_order), by = "Indiv")
  
  pop_levels <- unique(admixed$Pop[order(admixed$Pop)])
  
  admixed <- admixed %>%
    arrange(factor(Pop, levels = pop_levels), within_pop_order) %>%
    mutate(Indiv = factor(Indiv, levels = unique(Indiv)))
  
  #### Plot #################################################
  # population bar
  
  ## add a column called Locality_text that just has Locality as all the values (x-axis)
  admixed <- admixed %>% 
    mutate(Locality_text = "Locality")
  
  # admix data
  admix_adaptive_plot[[K]] <- ggplot(admixed, aes(x = Indiv, y = Proportion, fill = Cluster)) +
    geom_bar(width = 1, stat = "identity", position = "stack") +
    facet_grid( ~ Pop, scales = "free_x", space = "free_x") +
    labs(y = paste0("K",K), fill = "Cluster") +
    theme_minimal() +
    theme(axis.line = element_blank(), 
          title = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          plot.margin = unit(c(t = 0,r = 0,l = 0,b = 0), "cm"),
          axis.text.y = element_text(size = 12, color = 'black'),
          axis.text.x=element_blank(), 
          axis.ticks.y = element_line(),
          axis.title.y = element_text(hjust=1, angle=0, size = 18),
          axis.title.x = element_blank(), strip.text.x = element_blank()) +
    scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1),
                       expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = admix_colors)
  
}

admix_adaptive_plot[-length(admix_adaptive_plot)] <- lapply(admix_adaptive_plot[-length(admix_adaptive_plot)], 
                                                          function(p) p + theme(legend.position = "none"))

adaptive_admix_K24 <- Reduce('/',Filter(Negate(is.null), admix_adaptive_plot)) / pop_plot + 
  plot_layout(guides = "collect",
              heights = c(rep(1,length(Filter(Negate(is.null), admix_adaptive_plot))),0.25)) +
  plot_annotation(title = 'Fused: Adaptive Loci',
                  theme = theme(plot.title = element_text(hjust=0.1, size=20,
                                                          margin = margin(0,0,10,0))))
adaptive_admix_K24 


## PCAs #######################################################################

WILDCARD = '*fuse*'
COVFILES <- Sys.glob(file.path(here(),"results","pca",WILDCARD))
COVFILES <- COVFILES[!str_detect(COVFILES,pattern = 'filter')] # remove subset PCAs
COVFILES <- COVFILES[!str_detect(COVFILES,pattern = 'prune')] # remove pruned PCAs
COVFILES <- COVFILES[!str_detect(COVFILES,pattern = 'unfused_adap')] # remove unfused adaptive PCAs
  basename(COVFILES)

popFID <- left_join(bam_df, pop_df, by = "sampleID")

########### Read in PCA #################################
pca.vectors.list <- list()
varPC1 <- list(); varPC2 <- list()

for(i in 1:length(COVFILES)){
  
  # read in the covariance matrix
  cov <- as.matrix(read.table(COVFILES[i]))
  e <- eigen(cov)                            # calculate eigenvector values
  e_vectors <- as.data.frame(e$vectors) %>% mutate(FID = row_number())  # add FID to e_vectors
  e_per <- e$values/sum(e$values)            # percent explained by each component
  
  # determine the variance explained as a percent
  varPC1[[i]] <- (e$values[1]/sum(e$values))*100 #Variance explained by PC1
  varPC2[[i]] <- (e$values[2]/sum(e$values))*100 #Variance explained by PC2
  
  ##combine row names (population info) with the covariance matrix
  pca.vectors.list[[i]] = as_tibble(cbind(popFID, e_vectors))[,1:8]
  
}      

pca.vectors.list$Pop <- factor(pca.vectors.list$Pop, levels = color_df$Pop)

##### Plot theme #######################

COVFILES
# manually made to match order of covfiles
TITLES <- c("Fused: Adaptive Loci", "Fused: All Loci", 
            "Unfused: All Loci", "Unfused: Neutral Loci") 

theme_set(
  theme(legend.title=element_text(size = 15), legend.text=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(angle = 0, size = 10), axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"), panel.spacing = unit(0,"lines"),
        legend.position = "right", strip.text.y = element_text(angle = 0),
        title = element_text(size = 12)
  )
)

pc12 <- list()

for(j in seq(length(TITLES))) {  # subtract 1 because of factoring which adds up to 5
  
  pc12[[j]] <- as.data.frame(pca.vectors.list[[j]]) %>%
    mutate(Pop = fct_relevel(Pop, levels(pca.vectors.list$Pop))) %>%
    ggplot(aes(x=V1, y=V2, fill = Pop)) + 
    geom_point(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
    scale_fill_manual(name = "Population", values = mypalette) +
    ggtitle(TITLES[j]) +
    labs(x = paste0("PC1 (",round(varPC1[[j]], digits = 2),"%)"), 
         y= paste0("PC2 (",round(varPC2[[j]], digits = 2),"%)")) +
    theme(legend.position = 'none')
  
}


### ADMIX & PCA Plot ###########################################################

adaptive_pcas <- (pc12[[2]] / plot_spacer() / pc12[[1]]) + plot_layout(heights = c(3, 0.1, 3))
adaptive_pcas

adaptive_plots <- (adaptive_pcas | plot_spacer() | adaptive_admix_K24) + 
  plot_layout(guides = "collect", widths = c(3, 0.1, 6))
adaptive_plots

ggsave(plot = adaptive_plots, "./figures/figure4/boreogadus_adaptive_pca-admix_K4_20250707.jpeg",
       width = 15, height = 8)

neutral_pcas <- (pc12[[3]] / plot_spacer() / pc12[[4]]) + plot_layout(heights = c(3, 0.1, 3))

neutral_plots <- (neutral_pcas | plot_spacer() | neutral_admix_K23) + 
  plot_layout(guides = "collect", widths = c(3, 0.1, 6))
neutral_plots

ggsave(plot = neutral_plots, "./figures/figure4/boreogadus_neutral_pca-admix_K3_20250707.jpeg",
       width = 15, height = 8)

neutral_plots_K4 <- (neutral_pcas | plot_spacer() | neutral_admix_K24) + 
  plot_layout(guides = "collect", widths = c(3, 0.1, 6))
neutral_plots_K4

ggsave(plot = neutral_plots_K4, "./figures/figure4/boreogadus_neutral_pca-admix_K4_20250707.jpeg",
       width = 15, height = 8)

# I couldnt get a good final plot to work, so moving to Inkscape and Finalizing it there
# final_plots <- ((neutral_plots / plot_spacer() / adaptive_plots)) + 
#   plot_layout(guides = "collect", heights = c(5, 0.1, 5)) +
#   plot_annotation(tag_levels = c('A','B','C','','','D','E','F','','')) &
#   theme(plot.tag = element_text(face = 'bold', size = 22))
# final_plots

