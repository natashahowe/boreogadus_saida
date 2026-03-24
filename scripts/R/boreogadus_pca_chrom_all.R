


library(here)
library(tools)
library(tidyverse)
library(patchwork)

#### OPTIONS #########################################################
suffix = '_20241222'
WILDCARD = '*OZ*.1_2024*'
COVFILES <- Sys.glob(file.path(here(),"results","pca",WILDCARD))
basename(COVFILES)

## Filepaths ################################################################
BAMFILE <- "./data/R/boreogadus_filtered_downsampled_bams.txt"
METADATAFILE <- "./data/R/boreogadus_metadata.csv"
COLORFILE <- "./data/R/color_metadata_allpops_flip.txt"
chrom_df <- read.table('./data/R/boreogadus_chromosomes.txt', header = T)
filter_df <- as.data.frame(read.delim2("./data/boreogadus_unfused_chroms_prefilter.txt", header = F))

### Metadata ###############################################################
meta_df <- read.csv(METADATAFILE, header = T) %>%
  mutate(Population = ifelse(Population=='Southhampton Island', 'Southampton Island', Population),
         Pop = case_when(Population == 'SE Baffin' ~ 'SEBaffin',
                         TRUE ~ gsub(' .*$','',Population))) %>%
  filter(Pop != 'Unknown')

color_df <- read.delim2(COLORFILE, header = T, row.names = NULL, sep = "\t") %>%
  mutate(Pop = gsub('hh','h',Pop)) %>%
  left_join(meta_df %>% distinct(Population, Pop))

# call in bams
bam_df <- read.table(BAMFILE, header = F) %>%
  mutate(sampleID = basename(V1) %>% tools::file_path_sans_ext(.),
         sampleID = gsub("_sorted",'',sampleID) %>% gsub("_downsampled",'',.),
         sampleID = gsub("boreogadus_",'',sampleID)) %>%
  select(sampleID)

# call in some metadata
pop_df <- meta_df %>% # remove everything after space
  left_join(bam_df, ., by = "sampleID") 

pop_df <- cbind(pop_df, filter_df) %>% rename(Filter = V1)

########### read in PCA #################################
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
  #pca.vectors.list[[i]] = as_tibble(cbind(popFID, e_vectors))[,1:8]
  pca.vectors.list[[i]] = as_tibble(cbind(pop_df, e_vectors))[,1:8]
  
}      

#### PLOT THEME #######################

theme_set(
  theme(legend.title=element_blank(), legend.text=element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(angle = 0, size = 8), axis.title = element_text(size = 9),
        panel.background = element_rect(fill = "white"), panel.spacing = unit(0,"lines"),
        legend.position = "right", strip.text.y = element_text(angle = 0),
        title = element_text(size = 10)
  )
)

mypalette <- color_df$Color2
names(mypalette) <- color_df$Pop

pc12 <- list()

for(j in seq(length(pca.vectors.list))) {  # subtract 1 because of factoring which adds another value for some reason
  
  pca.vectors.i <- as.data.frame(pca.vectors.list[[j]])
  pca.vectors.i$Pop <- factor(pca.vectors.i$Pop, levels = color_df$Pop, ordered = T)
  
  file_name = gsub('boreogadus_','',basename(file_path_sans_ext(COVFILES[j])))
  chr_name = gsub(suffix,'',file_name)
  chr = chrom_df$chr[which(chrom_df$chrName == chr_name)]
  
  pc12[[j]] <- ggplot(data = pca.vectors.i, 
                      aes(x=V1, y=V2, fill = Pop)) + 
    geom_point(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
    scale_fill_manual(values = mypalette) +
    ggtitle(paste('Chr',chr)) +
    labs(x = paste0("PC1 (",round(varPC1[[j]], digits = 2),"%)"), 
         y= paste0("PC2 (",round(varPC2[[j]], digits = 2),"%)")) +
    scale_x_continuous(expand = expansion(0.05,0.05)) +
    scale_y_continuous(expand = expansion(0.05,0.05))
  
}

pcs12 <- Reduce('+',pc12) + plot_layout(guides = "collect",
                                        nrow = 3)
pcs12

ggsave(paste0("./figures/pca/20241222/boreogadus_chroms",suffix,"_allIndivs.jpg"),
       plot = pcs12, width = 16, height = 8, units = "in")

pcs12_tall <- Reduce('+',pc12) + plot_layout(guides = "collect",
                                        ncol = 4)

ggsave(paste0("./figures/pca/20241222/boreogadus_chroms",suffix,"_allIndivs_Tall.jpg"),
       plot = pcs12_tall, width = 10, height = 11, units = "in")

#### REMOVE FOUR TROUBLESOME INDIVIDUALS

pc12_filt <- list()

for(j in seq(length(pca.vectors.list))) {  # subtract 1 because of factoring which adds another value for some reason
  
  pca.vectors.i <- as.data.frame(pca.vectors.list[[j]]) %>%
    filter(Filter != 0)
  pca.vectors.i$Pop <- factor(pca.vectors.i$Pop, levels = color_df$Pop, ordered = T)
  
  file_name = gsub('boreogadus_','',basename(file_path_sans_ext(COVFILES[j])))
  chr_name = gsub(suffix,'',file_name)
  chr = chrom_df$chr[which(chrom_df$chrName == chr_name)]
  
  pc12_filt[[j]] <- ggplot(data = pca.vectors.i, 
                           aes(x=V1, y=V2, fill = Pop)) + 
    geom_point(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
    scale_fill_manual(values = mypalette) +
    ggtitle(paste('Chr',chr)) +
    labs(x = 'PC1', y= 'PC1') +
    scale_x_continuous(expand = expansion(0.02,0.02)) +
    scale_y_continuous(expand = expansion(0.02,0.02))
  
}

pca12_filt <- Reduce('+',pc12_filt) + plot_layout(guides = "collect")
pca12_filt

ggsave(paste0("./figures/pca/20241222/fusion/boreogadus_chroms",suffix,"_removeFour.jpg"),
       plot = pca12_filt, width = 12, height = 12, units = "in")

