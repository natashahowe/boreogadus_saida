############################################################################
## Genotype Likelihood Heatmap
# Adapted from Sara Schaal's R script for genotype likelihood data

############################################################################
### INSTALL PACKAGES & LOAD FUNCTIONS

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "data.table", 
                     "here", "viridis", "patchwork")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

#################################################################################
prefix = "boreogadus"

bglname="boreogadus_OZ177908.1_s26000000_e27500000" # beagle file name

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

# # create mypalette
mypalette <- as.vector(color_df$Color2) # turn colors into vector
  names(mypalette) <- color_df$Pop # attach pop name to palette color
  pop.factor.levels <- color_df$Pop

### Beagle and Cov #########################################

# Beagle file that you are interested in plotting the genotype matrix for
beagle_file <- read.table(paste0("./results/gls/",bglname,".beagle.gz"), header = F)

cov <- as.matrix(read.table(paste0("./results/pca/",bglname,".cov")))
  e <- eigen(cov)
  e_vectors <- as.data.frame(e$vectors) 

order_PC1 <- e_vectors %>%
  mutate(Ind = row_number()-1) %>%
  select(Ind, V1)

order_PC1 <- order_PC1[order(order_PC1$V1),] # this will be to check that the individuals are ordered according to pc1

###################################

# join those two dataframes
plotData <- pop_df %>%
  select(Pop)

## add an individual identifier that will match your genotype matrix labels
plotData$Ind <- paste0("Ind", 0:(nrow(plotData)-1))
head(plotData)

### Add PC1 data for gradient
plotData <- cbind(plotData, (e_vectors %>% select(V1)))
colnames(plotData) <- c("Pop", "Ind", "PC1")
head(plotData)

#################################################################################

plot.matrix <- matrix(ncol = ((ncol(beagle_file)-3)/3), nrow = nrow(beagle_file))
rownames(plot.matrix) <- beagle_file$V1
colnames(plot.matrix) <- paste0("Ind",0:(((ncol(beagle_file)-3)/3)-1))
for(i in 1:nrow(beagle_file)){
  for(j in 4:ncol(beagle_file)){
    if(j%%3 == 0){
      
      hom_major <- beagle_file[i,(j-2)]
      het <- beagle_file[i,(j-1)]
      hom_minor <- beagle_file[i,(j)]
      probs <- c(hom_major, het, hom_minor)
      
      if(probs[1] == probs[2] & probs[1] == probs[3]){ ## all equal make WHITE
        plot.matrix[i,((j/3)-1)] <- NA
      }else if(probs[1] == max(probs) & probs[1] == probs[3]){ ## SHOULD NEVER HAPPEN
        plot.matrix[i,((j/3)-1)] <- "NO"
      }else if(probs[1] == max(probs) & probs[1] == probs[2]){ ## homD==het --> RANGE 0.25:0.5 (actual range is 0.445:0.5) but this is a rare occurrence
        plot.matrix[i,((j/3)-1)] <- (1 + probs[1])/3
      }else if(probs[2] == max(probs) & probs[2] == probs[3]){## homR==het --> RANGE -0.25:-0.5 (actual range is -0.445:-0.5) but this is rare occurrence
        plot.matrix[i,((j/3)-1)] <- -1*(1 + probs[2])/4
      }else if(max(probs) == probs[1]){ ## homD --> RANGE 0.5:1
        plot.matrix[i,((j/3)-1)] <- (1 + probs[1])/2 
      }else if(max(probs) == probs[2] & probs[1] > probs[3]){ ## het, leaning homD --> RANGE 0:0.25
        plot.matrix[i,((j/3)-1)] <- (1 - probs[2])/2
      }else if(max(probs) == probs[2] & probs[1] < probs[3]){ ## het, leaning homR --> RANGE 0:-0.25
        plot.matrix[i,((j/3)-1)] <- (probs[2] - 1)/2
      }else if(max(probs) == probs[2] & probs[1] == probs[3]){ ## het, no lean --> VALUE 0
        plot.matrix[i,((j/3)-1)] <- 0 
      }else if(max(probs) == probs[3]){ ## homR --> RANGE -0.5:-1
        plot.matrix[i,((j/3)-1)] <- -1*(1 + probs[3])/2
      }
    }
  } 
}

## Take matrix and create a data.table for plotting
plot.table <- reshape2::melt(plot.matrix)

head(plot.table)

## Rename columns 
colnames(plot.table) <- c("locus", "Ind", "GenotypeScore")

plot.table <- plot.table %>%
  mutate(locus = gsub('.*.1_','',locus),
         Ind = as.character(Ind))
head(plot.table)

# ## Join your metadata with the new dataframe with genotype scores
# mypalette <- c("red3", "darkorange","cyan4")
# names(mypalette) <- c('Sex1a', 'Sex1b', 'Sex2')
# mypalette

plot.table.meta <- left_join(plot.table, plotData, by = "Ind") %>%
  mutate(Pop = factor(Pop, levels = pop.factor.levels))

GENOTITLE <- "Population"

geno_palette <- c("dodgerblue3","#21918c", "#5ec962", "#9be52a", "goldenrod1")

###################################################################

# #order for plotting
# plot.table.meta$Ind <- factor(plot.table.meta$Ind,
#                               levels=unique((plot.table.meta$Ind)[order(plot.table.meta$Pop)]))

#order based on PC1
plot.table.meta$Ind <- factor(plot.table.meta$Ind,
                              levels=unique((plot.table.meta$Ind)[order(plot.table.meta$Pop, -plot.table.meta$PC1)]))

## Make the colored bar that is for your population colors 
## add a column called Locality_text that just has Locality as all the values 
## this will serve as your X-axis and your population ID will serve as your y-axis
plot.table.meta <- plot.table.meta %>% 
  mutate(Locality_text = "Locality")
head(plot.table.meta)

plot.table.meta$locus <- as.character(plot.table.meta$locus)

############################################################################

# Plot your genotype heatmap 

geno_heatmap <- ggplot(plot.table.meta,aes(x=locus, y=Ind, fill=GenotypeScore)) + 
  geom_tile() + 
  theme_minimal()+ 
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(linewidth = 0.8),
        text = element_text(size = 11), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = "right", 
        legend.key.size = unit(0.6, 'cm'),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_blank())+
  scale_fill_gradientn(name = "Genotype Likelihood", 
                       colors = geno_palette,
                       guide = "colourbar", na.value = "gray90",
                       breaks = c(-1, -0.5, 0, 0.5, 1), 
                       labels = c("P(Alt) = 1", "P(Het) = P(Alt)", "P(Het) = 1", 
                                  "P(Ref) = P(Het)", "P(Ref) = 1"),
                       limits = c(-1,1)) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 5/1) 

# plot a bar that represents population ids for the individuals in rows of the genotype heatmap
pop_heatmap <- ggplot(plot.table.meta, aes(Locality_text,Ind, fill=Pop)) + 
  geom_tile() + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),  axis.line = element_blank(), 
        legend.position = "right", legend.title = element_text(size = 12),
        legend.text = element_text(size = 11), legend.key.size = unit(0.6, 'cm'),
        plot.margin = unit(c(0,0,0,0), "cm"))+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank())+
  scale_fill_manual(values = mypalette, name = GENOTITLE) +
  coord_fixed(ratio = 1/10) + 
  theme(axis.text.x=element_blank()) 

geno_plot <- (geno_heatmap + pop_heatmap) + plot_layout(guides = "collect") & theme(legend.position = "right")
#geno_plot
jpeg(paste0("./figures/genotype/",bglname,"_byPop_GL.jpeg"), width = 20, height = 5, res = 150, units = "in")
print(geno_plot)
dev.off()

#### END OF SCRIPT