# top code is for individual PCAs
# lower code is for local PCAs

library(tidyverse)
library(ggrepel)
library(here)
library(tools)
library(stringr)
here()

meta_df <- read.csv("./data/R/boreogadus_metadata.csv", header = T) %>%
  mutate(Population = ifelse(Population == 'Southhampton Island','Southampton Island',Population),
         Pop = case_when(Population == 'SE Baffin' ~ 'SEBaffin',
                         TRUE ~ gsub(' .*$','',Population))) %>%
  filter(Pop != 'Unknown')
meta_df

# read in color file
color_df <- read.delim("./data/R/color_metadata_allpops_flip.txt", header = T) %>%
  mutate(Pop = gsub('hh','h',Pop)) %>%
  left_join(meta_df %>% distinct(Population,Pop))  
color_df

#########################################################################

suffixes <- c("20241222", "noChukchi_20241222", "noChukchiCoronation_20241222", 
              "noChukchiIceland_20241222", "noIceland_20241222")

suffix <- "20241222"
#for(i in 1:5){
#  suffix = suffixes[i]

  # read in the covariance matrix
  cov <- read.table(paste0("./results/pca/boreogadus_wholegenome_",suffix,".cov"))
  
    # calculate eigenvector values
    e <- eigen(cov)
    e_vectors <- as.data.frame(e$vectors)
    e_per <- e$values/sum(e$values) # percent explained by each component
    
  # Get ID and pop info for each individual
  
    FID <- read.table(paste0("./data/bam/boreogadus_filtered_downsampled_bams.txt"),
                    sep = "\t", header = F) %>%
      mutate(temp = basename(file_path_sans_ext(V1)),
             temp = sub("_downsampled","",temp),
             temp = sub("boreogadus_","",temp),
             sampleID = sub("_sorted","",temp)) %>%
      left_join(meta_df, by = "sampleID") %>% 
      mutate(Population = ifelse(is.na(Population),"Iceland",Population)) %>%
      dplyr::select(sampleID, Population)
      
  #suffixes without pops
    if(i==2){FID <- FID %>% filter(Population != "Chukchi Sea")
    }else if(i==3){FID <- FID  %>% filter(Population != "Chukchi Sea",
                                          Population != "Coronation Gulf")
    }else if(i==4){FID <- FID  %>% filter(Population != "Chukchi Sea",
                                          Population != "Iceland")
    }else if(i==5){FID <- FID  %>% filter(Population != "Iceland")}
    
  #combine row names (population info) with the covariance matrix
  pca.vectors = as_tibble(cbind(FID, e_vectors)) 
  # determine the variance explained as a percent
  pca.eigenval.sum = sum(e$values) #sum of eigenvalues
  varPC1 <- format(round((e$values[1]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC1
  varPC2 <- format(round((e$values[2]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC2
  varPC3 <- format(round((e$values[3]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC3
  varPC4 <- format(round((e$values[4]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC4
  varPC5 <- format(round((e$values[5]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC4
  
  ###############################################################################
  
  # create mypalette
  mypalette <- as.vector(color_df$Color2)
  names(mypalette) <- color_df$Population # attach pop name to palette color
  pop.factor.levels <- color_df$Population
  mypalette
  
  pca.vectors$Population <- factor(pca.vectors$Population, levels = pop.factor.levels)
  levels(pca.vectors$Population)
  
  #levels(pca.vectors$Population)[levels(pca.vectors$Population)=='Southhampton'] <- 'Southampton'
  
  theme_set(
    theme( 
      legend.title=element_text(size = 14, hjust=0), 
      legend.text=element_text(size=11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(angle = 0, size = 10),
      axis.title = element_text(size = 14),
      panel.background = element_rect(fill = "white"), 
      panel.spacing = unit(0,"lines"),
      strip.text.y = element_text(angle = 0)
    )
  )
  
  pca1_2 <- ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) + 
    #geom_text_repel() +
    labs(x = paste0("PC1 (",varPC1,"% Variance)"), y= paste0("PC2 (",varPC2,"% Variance)"))
  pca1_2
  ggsave(filename = paste0("./figures/pca/boreogadus_wholegenome_pca1_2_",suffix,".jpeg"),
         plot = pca1_2, width = 9, height = 7, units = "in")
  
  ### ADDITIONAL PCAs ##################
  pca2_3 <- ggplot(data = pca.vectors, aes(x=V2, y=V3, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) +  
    #geom_text_repel() +
    labs(x = paste0("PC2 (",varPC2,"% Variance)"), y= paste0("PC3 (",varPC3,"% Variance)")) 
  pca2_3
  
  pca1_3 <- ggplot(data = pca.vectors, aes(x=V1, y=V3, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) +  
    #geom_text_repel() +
    labs(x = paste0("PC1 (",varPC1,"%)"), y= paste0("PC3 (",varPC3,"%)")) 
  pca1_3

  pca1_4 <- ggplot(data = pca.vectors, aes(x=V1, y=V4, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) +  
    labs(x = paste0("PC1 (",varPC1,"% Variance)"), y= paste0("PC4 (",varPC4,"% Variance)"))
  pca1_4
    
  pca2_4 <-  pca.vectors %>%
    filter(V4 < 0.4) %>%
    ggplot(., aes(x=V2, y=V4, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) +  
    labs(x = paste0("PC2 (",varPC2,"% Variance)"), y= paste0("PC4 (",varPC4,"% Variance)"))
  pca2_4
  
  pca2_4 <- pca.vectors %>%
    filter(V4 < 0.4) %>%
    ggplot(., aes(x=V2, y=V4, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) +
    labs(x = paste0("P2 (",varPC2,"% Variance)"), y= paste0("PC4 (",varPC4,"% Variance)"))
  pca2_4
  
  pca3_4 <-  pca.vectors %>%
    filter(V3 > -0.3, V4 < 0.4) %>%
    ggplot(., aes(x=V3, y=V4, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) +
    labs(x = paste0("P3 (",varPC3,"% Variance)"), y= paste0("PC4 (",varPC4,"% Variance)"))
  pca3_4
  
  pca2_5 <- pca.vectors %>%
    filter(V5 > -0.3) %>%
    ggplot(., aes(x=V2, y=V5, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) +
    labs(x = paste0("P2 (",varPC2,"% Variance)"), 
         y= paste0("PC5 (",varPC5,"% Variance)"))
  pca2_5
  
  ggsave(filename = paste0("./figures/pca/boreogadus_wholegenome_pca2_3_",suffix,".jpeg"), 
         plot = pca2_3, width = 10, height = 8, units = "in")
  ggsave(filename = paste0("./figures/pca/boreogadus_wholegenome_pca1_3_",suffix,".jpeg"),
         plot = pca1_3, width = 10, height = 8, units = "in")
  ggsave(filename = paste0("./figures/pca/boreogadus_wholegenome_pca1_4_",suffix,".jpeg"),
         plot = pca1_4, width = 10, height = 8, units = "in")
}


#### CHROM LEVEL PLOTS FOR ALL POPS #####

suffix <- "20241222"

for(i in 1:nrow(chrom_df)){
  
  chromName = chrom_df$chrName[i]
  chrom = chrom_df$chr[i]
  
  # read in the covariance matrix
  cov <- read.table(paste0("./results/pca/chr/boreogadus_",chromName,"_",suffix,".cov"))
  
  # calculate eigenvector values
  e <- eigen(cov)
  e_vectors <- as.data.frame(e$vectors)
  e_per <- e$values/sum(e$values) # percent explained by each component
  
  # Get ID and pop info for each individual
  
  FID <- read.table(paste0("./data/bam/boreogadus_filtered_downsampled_bams.txt"),
                    sep = "\t", header = F) %>%
    mutate(temp = basename(file_path_sans_ext(V1)),
           temp = sub("_downsampled","",temp),
           temp = sub("boreogadus_","",temp),
           sampleID = sub("_sorted","",temp)) %>%
    left_join(meta_df, by = "sampleID") %>% 
    mutate(Population = ifelse(is.na(Population),"Iceland",Population)) %>%
    dplyr::select(sampleID, Population)
  
  #combine row names (population info) with the covariance matrix
  pca.vectors = as_tibble(cbind(FID, e_vectors))
  # determine the variance explained as a percent
  pca.eigenval.sum = sum(e$values) #sum of eigenvalues
  varPC1 <- format(round((e$values[1]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC1
  varPC2 <- format(round((e$values[2]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC2
  varPC3 <- format(round((e$values[3]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC3

  ###############################################################################
  
  # create mypalette
  mypalette <- as.vector(color_df$Color2)
  names(mypalette) <- color_df$Population # attach pop name to palette color
  pop.factor.levels <- color_df$Population
  mypalette
  
  pca.vectors$Population <- factor(pca.vectors$Population, levels = pop.factor.levels)
  levels(pca.vectors$Population)
  
  theme_set(
    theme( 
      legend.title=element_blank(), 
      legend.text=element_text(size=11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(angle = 0, size = 10),
      axis.title = element_text(size = 14),
      panel.background = element_rect(fill = "white"), 
      panel.spacing = unit(0,"lines"),
      strip.text.y = element_text(angle = 0)
    )
  )
  
  pca1_2 <- ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) + 
    ggtitle(paste0("Boreogadus - Chr: ",chrom)) +
    labs(x = paste0("PC1 (",varPC1,"%)"), 
         y= paste0("PC2 (",varPC2,"%)"))
  pca1_2
  ggsave(filename = paste0("./figures/pca/chr/pc12/boreogadus_chr",chrom,"_pca1_2_",suffix,".jpeg"),
         plot = pca1_2, width = 10, height = 8, units = "in")
  
  ### ADDITIONAL PCAs ##################
  pca2_3 <- ggplot(data = pca.vectors, aes(x=V2, y=V3, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) +  
    ggtitle(paste0("Boreogadus - Chr: ",chrom)) +
    labs(x = paste0("PC2 (",varPC2,"%)"), 
         y= paste0("PC3 (",varPC3,"%)")) 
  pca2_3
  ggsave(filename = paste0("./figures/pca/chr/pc23/boreogadus_chr",chrom,"_pca2_3_",suffix,".jpeg"),
         plot = pca2_3, width = 10, height = 8, units = "in")
  
}


#### LOCAL PCAS #####

pca_list <- noquote(list.files(path = "./results/pca/peaks/",
                               pattern = "boreogadus_chr",
                               full.names = F))
  pca_list <- pca_list[grep(suffix, pca_list)]

for(i in 1:length(pca_list)){
  
  PCAFILE <- file_path_sans_ext(pca_list[[i]])
  # read in the covariance matrix
  cov <- read.table(paste0("./results/pca/peaks/",PCAFILE,".cov"))
  
  # calculate eigenvector values
  e <- eigen(cov)
  e_values <- e$values
  e_vectors <- as.data.frame(e$vectors)
  e_per <- e$values/sum(e$values) # percent explained by each component
  
  # Get ID and pop info for each individual
  
  FID <- read.table(paste0("./data/bam/boreogadus_filtered_downsampled_bams.txt"),
                    sep = "\t", header = F) %>%
    mutate(temp = basename(file_path_sans_ext(V1)),
           temp = sub("_downsampled","",temp),
           temp = sub("boreogadus_","",temp),
           sampleID = sub("_sorted","",temp)) %>%
    left_join(meta_df, by = "sampleID") %>% 
    mutate(Population = ifelse(is.na(Population),"Iceland",Population)) %>%
    dplyr::select(sampleID, Population)
  
  #combine row names (population info) with the covariance matrix
  pca.vectors = as_tibble(cbind(FID, e_vectors))
  # determine the variance explained as a percent
  pca.eigenval.sum = sum(e$values) #sum of eigenvalues
  varPC1 <- format(round((e$values[1]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC1
  varPC2 <- format(round((e$values[2]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC2
  
  ###############################################################################
  
  # create mypalette
  mypalette <- as.vector(color_df$Color2)
  names(mypalette) <- color_df$Population # attach pop name to palette color
  pop.factor.levels <- color_df$Population
  mypalette
  
  pca.vectors$Population <- factor(pca.vectors$Population, levels = pop.factor.levels)
  levels(pca.vectors$Population)
  
  theme_set(
    theme( 
      legend.title=element_blank(), 
      legend.text=element_text(size=11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(angle = 0, size = 10),
      axis.title = element_text(size = 14),
      panel.background = element_rect(fill = "white"), 
      panel.spacing = unit(0,"lines"),
      strip.text.y = element_text(angle = 0)
    )
  )
  
  pca1_2 <- ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Population)) + 
    geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
    scale_fill_manual(values = mypalette) + 
    ggtitle(paste0("Boreogadus Populations: ",PCAFILE)) +
    labs(x = paste0("PC1 (",varPC1,"% Variance)"), y= paste0("PC2 (",varPC2,"% Variance)"))
  pca1_2
  ggsave(filename = paste0("./figures/pca/",suffix,"/peaks/",PCAFILE,"_pca1_2.jpeg"),
         plot = pca1_2, width = 10, height = 8, units = "in")
}
