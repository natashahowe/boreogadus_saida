

library(here)
library(tidyverse)

DATE = "20241222"
markertypes = c('neutral_1SD','adaptive_2SD')
#alphas = c(0.05,0.001)
alpha = 0.05

METADATAFILE <- "./data/R/boreogadus_metadata.csv"
BAMFILE <- "./data/R/boreogadus_filtered_downsampled_bams.txt"

color_df <- read.delim2(paste0(here(),"./data/R/color_metadata_allpops_flip.txt"), 
                        header = T, row.names = NULL, sep = "\t")

#### PLOT THEME 

theme_set(
  theme(legend.title=element_blank(), legend.text=element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(angle = 0, size = 10), axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"), panel.spacing = unit(0,"lines"),
        legend.position = "right", strip.text.y = element_text(angle = 0)
  )
)

########### read in PCA #################################
for(i in 1:length(markertypes)){
  markertype = markertypes[i]
  #for(j in 1:length(alphas)){
    
    #alpha = alphas[j]
    
    # read in the covariance matrix
    cov <- as.matrix(read.table(paste0("./results/pca/boreogadus_all_",markertype,"_alpha",alpha,".cov")))
      e <- eigen(cov)                            # calculate eigenvector values
      e_vectors <- as.data.frame(e$vectors) %>%
        mutate(FID = row_number())               # add FID to e_vectors
      e_per <- e$values/sum(e$values)            # percent explained by each component
      
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
    
    popFID <- left_join(bam_df, pop_df, by = "sampleID")
    
    ##combine row names (population info) with the covariance matrix
    pca.vectors = as_tibble(cbind(popFID, e_vectors))
    
    # determine the variance explained as a percent
    pca.eigenval.sum = sum(e$values) #sum of eigenvalues
      varPC1 <- (e$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
      varPC2 <- (e$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
      varPC3 <- (e$values[3]/pca.eigenval.sum)*100 #Variance explained by PC2
      
    #################
    # PLOTTING
    
    pca.vectors$Pop <- factor(pca.vectors$Pop, levels = color_df$Pop)
      mypalette <- color_df$Color2
      names(mypalette) <- levels(pca.vectors$Pop)

    pc12 <- ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Pop), color = 'gray20') + 
      geom_point(alpha = 0.7, size = 2.5, pch = 21) +
      scale_fill_manual(values = mypalette) +
      ggtitle(paste0("Boreogadus ",markertype," Markers (alpha=",alpha,")")) +
      labs(x = paste0("PC1 (",round(varPC1, digits = 2),"%)"), 
           y= paste0("PC2 (",round(varPC2, digits = 2),"%)")) 
    pc12
    ggsave(filename = paste0("./figures/pca/",DATE,"/2SD/boreogadus_all_",markertype,"_alpha",alpha,"_pc1v2.jpeg"), 
           plot = pc12, width = 10, height = 8, units = "in")
    
    pc13 <- ggplot(data = pca.vectors, aes(x=V1, y=V3, fill = Pop), color = 'gray20') + 
      geom_point(alpha = 0.7, size = 2.5, pch = 21) +
      scale_fill_manual(values = mypalette) +
      ggtitle(paste0("Boreogadus ",markertype," Markers (alpha=",alpha,")")) +
      labs(x = paste0("PC1 (",round(varPC1, digits = 2),"%)"), 
           y= paste0("PC3 (",round(varPC3, digits = 2),"%)")) 
    pc13
    ggsave(filename = paste0("./figures/pca/",DATE,"/2SD/boreogadus_all_",markertype,"_alpha",alpha,"_pc1v3.jpeg"), 
           plot = pc13, width = 10, height = 8, units = "in")
    
  #}
}

