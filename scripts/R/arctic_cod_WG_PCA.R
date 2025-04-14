# top code is for individual PCAs
# lower code is for local PCAs

library(tidyverse)
library(ggrepel)

# choose plot naming convention

NAME = "10-03-27"
NAME = "166inds_10-03-27"
NAME = "166inds_10-26-23"

#########################################################################
# read in the covariance matrix
cov <- read.table(paste0("./results/pca/cod_wholegenome_",NAME,".cov"))

fid_num <- as.numeric(nrow(cov))

# calculate eigenvector values
e <- eigen(cov)
e_values <- e$values
e_vectors <- as.data.frame(e$vectors)

# percent explained by each component
e_per <- e$values/sum(e$values)

# Get ID and pop info for each individual
FID <- read.table("./data/R/COD_FID_166_pop_region.txt",
                  sep = "\t", header = T)

FID <- FID %>%
  filter(Pop != "Unknown",
         Pop != "Chukchi")

#combine row names (population info) with the covariance matrix
pca.vectors = as_tibble(cbind(FID, e_vectors))

# determine the variance explained as a percent
pca.eigenval.sum = sum(e$values) #sum of eigenvalues
varPC1 <- format(round((e$values[1]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC1
varPC2 <- format(round((e$values[2]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC2
varPC3 <- format(round((e$values[3]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC3
varPC4 <- format(round((e$values[4]/pca.eigenval.sum*100), 2), nsmall = 2) #Variance explained by PC4

###############################################################################

# read in color file
color_df <- read.delim("./data/R/color_metadata_allpops_flip.txt", 
                       header = TRUE)

# create mypalette
mypalette <- as.vector(color_df$Color2) # turn colors into vector
names(mypalette) <- color_df$Pop # attach pop name to palette color
pop.factor.levels <- color_df$Pop

pca.vectors$Pop <- factor(pca.vectors$Pop, levels = pop.factor.levels)

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

pca1_2 <- ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Pop, label = FID)) + 
  geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
  scale_fill_manual(values = mypalette) + 
  #geom_text_repel() +
  labs(x = paste0("PC1 (",varPC1,"% Variance)"), y= paste0("PC2 (",varPC2,"% Variance)"))
pca1_2
ggsave(filename = paste0("./figures/pca/wholegenome_pca1_2_cod_",NAME,"_noChukchi.jpeg"), 
       plot = pca1_2, width = 10, height = 8, units = "in")

### ADDITIONAL PCAs ##################
pca2_3 <- ggplot(data = pca.vectors, aes(x=V2, y=V3, fill = Pop, label = FID)) + 
  geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
  scale_fill_manual(values = mypalette) +  
  #geom_text_repel() +
  labs(x = paste0("PC2 (",varPC2,"% Variance)"), y= paste0("PC3 (",varPC3,"% Variance)")) 
pca2_3

pca1_3 <- ggplot(data = pca.vectors, aes(x=V1, y=V3, fill = Pop, label = FID)) + 
  geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
  scale_fill_manual(values = mypalette) +  
  #geom_text_repel() +
  labs(x = paste0("PC1 (",varPC1,"%)"), y= paste0("PC3 (",varPC3,"%)")) 
pca1_3

pca1_4 <- ggplot(data = pca.vectors, aes(x=V1, y=V4, fill = Pop)) + 
  geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
  scale_fill_manual(values = mypalette) +  
  labs(x = paste0("PC1 (",varPC1,"% Variance)"), y= paste0("PC4 (",varPC4,"% Variance)"))
pca1_4

pca2_4 <- ggplot(data = pca.vectors, aes(x=V2, y=V4, fill = Pop)) + 
  geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
  scale_fill_manual(values = mypalette) +
  labs(x = paste0("P2 (",varPC2,"% Variance)"), y= paste0("PC4 (",varPC4,"% Variance)"))
pca2_4

pca3_4 <- ggplot(data = pca.vectors, aes(x=V3, y=V4, fill = Pop)) + 
  geom_point(alpha = 0.85, pch = 21, color = "gray20", size = 3) +
  scale_fill_manual(values = mypalette) +
  labs(x = paste0("P3 (",varPC3,"% Variance)"), y= paste0("PC4 (",varPC4,"% Variance)"))
pca3_4

ggsave(filename = paste0(DIRECT,"/pca/pca_plots/pca2_3_cod_",NAME,"_newcol.jpeg"), 
       plot = pca2_3, width = 10, height = 8, units = "in")
ggsave(filename = paste0(DIRECT,"/pca/pca_plots/pca1_3_cod_",NAME,"_newcol.jpeg"), 
       plot = pca1_3, width = 10, height = 8, units = "in")
ggsave(filename = paste0(DIRECT,"/pca/pca_plots/pca1_4_cod_",NAME,"_newcol.jpeg"), 
       plot = pca1_4, width = 10, height = 8, units = "in")
