# Chroms 6-18 - Neutral

library(here)
library(tools)
library(tidyverse)
library(patchwork)

#### OPTIONS #########################################################
WILDCARD = '*prun*1MB*'
#WILDCARD = '*prun*filt*1MB*'
COVFILES <- Sys.glob(file.path(here(),"results","pca",WILDCARD))
  COVFILES <- COVFILES[!str_detect(COVFILES,pattern = 'OZ')]
  #COVFILES <- COVFILES[!str_detect(COVFILES,pattern = 'filter')]
  basename(COVFILES)

METADATAFILE <- "./data/R/boreogadus_metadata.csv"
BAMFILE <- "./data/R/boreogadus_filtered_downsampled_bams.txt"

color_df <- read.delim2(paste0(here(),"./data/R/color_metadata_allpops_flip.txt"), 
                        header = T, row.names = NULL, sep = "\t")

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
                         TRUE ~ gsub(' .*$','',Population))) # remove everything after space

filter_df <- as.data.frame(read.delim2("./data/boreogadus_unfused_chroms_prefilter.txt", header = F))

popFID <- left_join(bam_df, pop_df, by = "sampleID") %>%
  cbind(filter_df) %>% rename(Filter = V1)

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
  varPC3[[i]] <- (e$values[3]/sum(e$values))*100 #Variance explained by PC3
  
  if(str_detect(COVFILES[i],pattern = 'filter')==F){
  ##combine row names (population info) with the covariance matrix
    pca.vectors.list[[i]] = as_tibble(cbind(popFID, e_vectors))[,1:8]
  }else{
    pca.vectors.list[[i]] = as_tibble(cbind(filter(popFID, Filter == 1), e_vectors))[,1:8]
  }
}      

#### PLOT THEME #######################

theme_set(
  theme(legend.title=element_blank(), legend.text=element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(angle = 0, size = 10), axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"), panel.spacing = unit(0,"lines"),
        legend.position = "right", strip.text.y = element_text(angle = 0)
  )
)

mypalette <- color_df$Color2
  names(mypalette) <- color_df$Pop

pc12 <- list(); pc13 <- list()
  
for(j in seq(length(pca.vectors.list))) {  # subtract 1 because of factoring which adds another value for some reason
  
  pca.vectors.i <- as.data.frame(pca.vectors.list[[j]])
  pca.vectors.i$Pop <- factor(pca.vectors.i$Pop, levels = color_df$Pop, ordered = T)
  
  pc12[[j]] <- ggplot(data = pca.vectors.i, 
                      aes(x=V1, y=V2, fill = Pop)) + 
    geom_point(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
    scale_fill_manual(values = mypalette) +
    ggtitle(gsub('boreogadus_unfused_neutral_','',basename(file_path_sans_ext(COVFILES[j])))) +
    labs(x = paste0("PC1 (",round(varPC1[[j]], digits = 2),"%)"), 
         y= paste0("PC2 (",round(varPC2[[j]], digits = 2),"%)")) 

}

pcas <- Reduce('+',pc12) / Reduce('+',pc13) + 
  plot_layout(guides = "collect")
pcas

pcs12 <- Reduce('+',pc12) + plot_layout(guides = "collect")
ggsave("./figures/pca/20241222/fusion/boreogadus_unfused_chroms_neutral_1MB.jpg",
       plot = pcs12, width = 12, height = 6, units = "in")

# remove Iceland

pca.vectors.2 <- as.data.frame(pca.vectors.list[[2]]) %>%
  filter(Pop != 'Iceland')
pca.vectors.2$Pop <- factor(pca.vectors.2$Pop, levels = color_df$Pop, ordered = T)

pcsNI <- ggplot(data = pca.vectors.2, 
                    aes(x=V1, y=V2, fill = Pop)) + 
  geom_point(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
  scale_fill_manual(values = mypalette) +
  ggtitle(gsub('boreogadus_unfused_neutral_','',basename(file_path_sans_ext(COVFILES[2])))) +
  labs(x = "PC1", y= "PC2") 

ggsave("./figures/pca/20241222/fusion/boreogadus_unfused_chroms_neutral_1MB_filterNoIceland.jpg",
       plot = pcsNI, width = 12, height = 6, units = "in")


#### before filtering individuals ##################################################

filterIndivs <- pca.vectors.list[[1]] %>%
  mutate(Filter = ifelse(V1 > 0.2, 0, 1))

write.table(filterIndivs[,"Filter"], "./data/boreogadus_unfused_chroms_prefilter.txt", 
            quote = F, col.names = F, row.names = F)
