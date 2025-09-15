

library(here)
library(tools)
library(tidyverse)
library(patchwork)

#### OPTIONS #########################################################
WILDCARD = '*fuse*'
COVFILES <- Sys.glob(file.path(here(),"results","pca",WILDCARD))
  COVFILES <- COVFILES[!str_detect(COVFILES,pattern = 'filter')]
  COVFILES <- COVFILES[!str_detect(COVFILES,pattern = 'prune')]
basename(COVFILES)

METADATAFILE <- "./data/R/boreogadus_metadata.csv"
BAMFILE <- "./data/R/boreogadus_filtered_downsampled_bams.txt"

color_df <- read.delim2(paste0(here(),"./data/R/color_metadata_allpops_flip.txt"), 
                        header = T, row.names = NULL, sep = "\t") %>%
  mutate(Pop = gsub('hh','h',Pop))

# call in bams
bam_df <- read.table(BAMFILE, header = F) %>%
  mutate(sampleID = basename(V1),
         sampleID = gsub("_sorted.bam",'',sampleID),
         sampleID = gsub("_sorted_downsampled.bam",'',sampleID),
         sampleID = gsub("boreogadus_",'',sampleID)) %>%
  select(sampleID)

# call in some metadata
pop_df <- read.csv(METADATAFILE, header = T) %>%
  mutate(Pop = case_when(Population == 'SE Baffin' ~ 'SEBaffin',
                         TRUE ~ gsub(' .*$','',Population)),  # remove everything after space
         Pop = gsub('hh','h',Pop))
  
popFID <- left_join(bam_df, pop_df, by = "sampleID")
unique(pop_df$Pop)

########### read in PCA #################################
pca.vectors.list <- list()
varPC1 <- list(); varPC2 <- list(); varPC3 <- list()

for(i in 1:length(COVFILES)){

  # read in the covariance matrix
  cov <- as.matrix(read.table(COVFILES[i]))
    e <- eigen(cov)                            # calculate eigenvector values
    e_vectors <- as.data.frame(e$vectors) %>% mutate(FID = row_number())  # add FID to e_vectors
    e_per <- e$values/sum(e$values)            # percent explained by each component
   
  # determine the variance explained as a percent
  varPC1[[i]] <- (e$values[1]/sum(e$values))*100 #Variance explained by PC1
  varPC2[[i]] <- (e$values[2]/sum(e$values))*100 #Variance explained by PC2
  varPC3[[i]] <- (e$values[3]/sum(e$values))*100 #Variance explained by PC2
  
  ##combine row names (population info) with the covariance matrix
  pca.vectors.list[[i]] = as_tibble(cbind(popFID, e_vectors))[,1:8]

}      

pca.vectors.list$Pop <- factor(pca.vectors.list$Pop, levels = color_df$Pop)
  mypalette <- color_df$Color2
  names(mypalette) <- levels(pca.vectors.list$Pop)

#### PLOT THEME #######################

COVFILES
# manually made to match order of covfiles
TITLES <- c("Fused Chromosomes: Adaptive Loci", "Fused Chromosomes: All Loci", 
            "Unfused Chromosomes: Adaptive Loci",
            "Unfused Chromosomes: All Loci", "Unfused Chromosomes: Neutral Loci") 
  
  
theme_set(
  theme(legend.title=element_text(size = 15), legend.text=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(angle = 0, size = 10), axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"), panel.spacing = unit(0,"lines"),
        legend.position = "right", strip.text.y = element_text(angle = 0)
  )
)

pc12 <- list(); pc13 <- list()

for(j in seq(length(COVFILES))) {  

  pc12[[j]] <- as.data.frame(pca.vectors.list[[j]]) %>%
    mutate(Pop = fct_relevel(Pop, levels(pca.vectors.list$Pop))) %>%
    ggplot(aes(x=V1, y=V2, fill = Pop)) + 
    geom_point(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
    scale_fill_manual(name = "Population", values = mypalette) +
    ggtitle(TITLES[j]) +
    labs(x = paste0("PC1 (",round(varPC1[[j]], digits = 2),"%)"), 
         y= paste0("PC2 (",round(varPC2[[j]], digits = 2),"%)")) 
  
}

# combine plots
pc12s <- (pc12[[2]] + pc12[[1]]) / 
  plot_spacer() / 
  (pc12[[4]] + pc12[[5]]) + 
  plot_layout(guides = "collect", 
              heights = c(2, 0.1, 2)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 22))
pc12s

ggsave("./figures/pca/20241222/fusion/boreogadus_fusion_chroms_adaptive-neutral_label.jpg",
       plot = pc12s, width = 12, height = 10, units = "in")

pc12[[3]]


ggsave("./figures/pca/20241222/fusion/boreogadus_unfused_chroms_adaptive.jpg",
       plot = pc12[[3]], width = 8, height = 6, units = "in")

### can add this to the j loop if want pca 1 v 3 #################
pc13[[j]] <- ggplot(data = as.data.frame(pca.vectors.list[[j]]), aes(x=V1, y=V3, fill = Pop)) + 
  geom_point(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
  scale_fill_manual(name = "Population", values = mypalette) +
  ggtitle(gsub('boreogadus_','',basename(file_path_sans_ext(COVFILES[j])))) +
  labs(x = paste0("PC1 (",round(varPC1[[j]], digits = 2),"%)"), 
       y= paste0("PC3 (",round(varPC3[[j]], digits = 2),"%)")) 

pc13s <- pc13[[2]] + pc13[[1]] + pc13[[4]] + pc13[[5]] + plot_layout(guides = "collect")
