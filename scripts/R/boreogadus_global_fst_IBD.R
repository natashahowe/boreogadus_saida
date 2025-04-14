######################################
### INSTALL PACKAGES & LOAD FUNCTIONS

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", 
                     "tools", "here", "reshape2")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

# CHOOSE OPTIONS
# for folder directory to angsd output and plot output
# for this code to work, will also need to match date in file names
# for example: Broughton.Frobisher_10-3-23.fst.SNP.txt
DATE="20241222" 

#################################################################################
# Read in file for order of plotting from West to East (kind of)
color_df <- read.delim("./data/R/color_metadata_allpops_flip.txt", header = T)

# create mypalette
mypalette <- as.vector(color_df$Color) # turn colors into vector
  names(mypalette) <- color_df$Pop # attach pop name to palette color
  mypalette

# read in marinedistances matrix (previously calculated & in km)  
  # distances were calculated with marinedistances, except Iceland
mat.geoDist <- as.matrix(read.csv("./results/ibd/geodistances_matrix_20250224.csv",
                         row.names = 1, header = T))
# all distances calculated in marinedistances, including Iceland
mat.geoDist <- as.matrix(read.csv("./results/ibd/boreogadus_noLabrador_marinedistances_res-h2.csv",
                                  row.names = 1, header = T))  

#################################################################################
# files for combining
  GLOBALNAMES <- Sys.glob(file.path(here(),"results","fst",DATE,"global","*20241222*")) 
  GLOBALNAMES <- GLOBALNAMES[-grep("_subset",GLOBALNAMES)] # not going to plot the subset pops in this
  GLOBALNAMES <- GLOBALNAMES[-grep("Labrador",GLOBALNAMES)]# too few individuals in labrador so not going to plot this round
  basename(GLOBALNAMES)
  global_list <- as.list(GLOBALNAMES)
  
  ################################################################################################
  # Read in pairwise comparison data files
  global_fst <- global_list %>%
    set_names(nm = GLOBALNAMES) %>%
    map_dfr(
      ~ read_delim(.x, col_types = cols(), 
                   col_names = c("unweighted", "weighted_fst"), delim = "\t"), .id = "comparison"
    )
  head(global_fst)
  
  # Strip the original file name such that only the population that FOCALPOP is 
  # being compared to is retained (for plot labeling)
  global_fst <- global_fst %>%
    mutate(temp = basename(tools::file_path_sans_ext(comparison)),
           temp = sub("boreogadus_","",temp),
           temp = sub("_subset","",temp),
           temp = sub("_downsampled","",temp),
           temp = sub(paste0("_",DATE,".sfs.global"),"",temp),
           POP1 = gsub("-.*","",temp),
           POP2 = gsub(".*-","",temp)) %>%
    select(c(POP1, POP2, weighted_fst))
  
  head(global_fst)
  unique(global_fst$POP2)
  
  # Make the POP2 column in the data frame into a factor with a specific order
  global_fst$POP1 <- factor(global_fst$POP1, levels = color_df$Pop)
  global_fst$POP2 <- factor(global_fst$POP2, levels = color_df$Pop)
  
  global_fst_rev <- global_fst %>%
    rename(POP2 = POP1, POP1 = POP2) %>% select(POP1,everything())

  mat.fst <- acast(rbind(global_fst,global_fst_rev), POP1~POP2, value.var = "weighted_fst")  

## Specify some functions to retrieve upper part of matrix
# Get lower triangle of the correlation matrix

get_lower_tri <- function(Fstmat){
  Fstmat[upper.tri(Fstmat)] <- NA
  return(Fstmat)
}

## subset the matrix
lower_tri <- get_lower_tri(mat.fst)
#View(lower_tri)

##Use the package reshape to melt the matrix into a df again:
final_df <- melt(lower_tri, value.name = "weighted_fst") %>%
  filter(!is.na(weighted_fst)) %>%
  mutate(weighted_fst = round(weighted_fst, digits = 4))

# order populations by the color_df dataframe order
final_df$POP1 <- factor(as.factor(final_df$Var1), levels = color_df$Pop)
final_df$POP2 <- factor(final_df$Var2, levels = color_df$Pop)

# Make a heatmap and visualize the FST values
heatmap_plot <- ggplot(data = final_df, aes(POP1, POP2, fill = weighted_fst)) +
  geom_raster() +
  geom_text(aes(label = weighted_fst), size = 3) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, vjust = 1, size = 14, hjust = 1),
        axis.text.y = element_text(color = "black", angle = 0, vjust = 0.5, size = 14, hjust = 1),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 15)) +
  ylab("Population A") +
  xlab("Population B") +
  labs(fill = expression(italic(F[ST]))) +
  coord_fixed() +
  theme(legend.position = c(0.2, .7))
heatmap_plot

# save plot to file
ggsave(paste0("./figures/fst/",DATE,"/global_Fst_heatmap2.jpeg"), 
       heatmap_plot, width= 7, height = 6)

################################################################################
### IBD ###

## convert matrices to a dataframe
df.fst <- melt(mat.fst) %>% filter(Var1 != Var2)
colnames(df.fst) <- c("POP1", "POP2", "weightedFST")

df.geoDist <- melt(mat.geoDist) %>% filter(Var1 != Var2)
colnames(df.geoDist) <- c("POP1", "POP2", "geoDist")

df.ibd <- full_join(df.fst, df.geoDist, by = c("POP1","POP2"))

# calculate genetic distance from FST
df.ibd$geneDist <- as.numeric(df.ibd$weightedFST)/(1-as.numeric(df.ibd$weightedFST))
df.ibd$POP1 <- factor(df.ibd$POP1, levels = color_df$Pop)

View(df.ibd)

### PLOTTING ###
# Set the ggplot theme
theme_set(
  theme( 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 70, size = 11, vjust = 0.5, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(angle = 0, size = 11, color = "black"),
    axis.title.y = element_text(size = 18, color = "black", angle = 90),
    title = element_text(size = 15, color = "black"),
    panel.background = element_rect(fill = "white"), 
    panel.spacing = unit(0,"lines"),
    strip.text.x = element_text(angle = 0, color = "black", size = 14),
    strip.text.y = element_text(angle = 0, color = "black", size = 13),
    strip.background = element_rect(fill = "gray90"), 
    axis.line = element_line(colour = "black")
  )
)

#### PLOT FACET GRID FOR EACH POPULATION ####
panelIBD <- ggplot(df.ibd, aes(x = geoDist, y = geneDist)) +
  geom_smooth(method = 'lm', se = F, col = "gray", alpha = 0.3) + 
  geom_point(aes(color = POP2), size = 3) + 
  facet_wrap(~POP1) +
  scale_color_manual(name = "Population", values = mypalette) + 
  stat_cor(aes(label = after_stat(r.label)), label.x = 0.1, label.y = 0.04)+
  stat_cor(aes(label = after_stat(rr.label)), label.x = 0.1, label.y = 0.048)+
  labs(title = "IBD in Arctic Cod",
       x = "Marine Distance (km)",
       y = expression("F"[ST]*"/(1-F"[ST]*")")) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 14), strip.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),legend.title = element_text(size = 14))
panelIBD
ggsave("./figures/ibd/boreogadus_allPops_ibd_facet_colorByPop_r_r2.jpeg",
       panelIBD, width = 10, height = 8, dpi = 300)

panelIBD <- ggplot(df.ibd, aes(x = geoDist, y = geneDist)) +
  geom_smooth(method = 'lm', se = F, col = "gray", alpha = 0.3) + 
  geom_point(aes(color = POP1), size = 3) + 
  facet_wrap(~POP1) +
  scale_color_manual(name = "Population", values = mypalette) + 
  #stat_cor(aes(label = after_stat(r.label)), label.x = 0.1, label.y = 0.046)+
  #stat_cor(aes(label = after_stat(rr.label)), label.x = 0.1, label.y = 0.054)+
  labs(title = "IBD in Arctic Cod",
       x = "Marine Distance (km)",
       y = expression("F"[ST]*"/(1-F"[ST]*")")) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 14), strip.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),legend.title = element_text(size = 14))
panelIBD
ggsave("./figures/ibd/boreogadus_allPops_ibd_facet_colorByPanel.jpeg",
       panelIBD, width = 10, height = 8, dpi = 300)

##### PLOT TOTAL IBD FOR ALL POP COMPARISONS #########

totalIBD <- ggplot(df.ibd, aes(x = geoDist, y = geneDist)) +
  geom_smooth(method = 'lm', se = F, col = "gray40", alpha = 0.3) + 
  geom_point(size = 4.5, pch = 21, stroke = 2, aes(fill = POP1, color = POP2)) +
  #stat_cor(aes(label = after_stat(r.label)), label.x = 200, label.y = 0.052)+
  #stat_cor(aes(label = after_stat(rr.label)), label.x = 200, label.y = 0.054)+
  scale_color_manual(name = "Population", values = mypalette) + 
  scale_fill_manual(name = "Population", values = mypalette) + 
  labs(title = "IBD in Arctic Cod",
       x = "Marine Distance (km)",
       y = expression("F"[ST]*"/(1-F"[ST]*")")) +
  theme_bw() + 
  theme(panel.border = element_blank(), axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 12),legend.title = element_text(size = 14))
totalIBD
ggsave("./figures/ibd/boreogadus_allPops_ibd_colorByPop.jpeg",
       totalIBD, width = 8, height = 6, dpi = 300)

# all points one color
totalIBD2 <- ggplot(df.ibd, aes(x = geoDist, y = geneDist)) +
  geom_smooth(method = 'lm', se = F, col = "gray60", alpha = 0.3) + 
  geom_point(size = 4, pch = 21, fill = "deepskyblue3", color = "black") +
  stat_cor(aes(label = after_stat(r.label)), label.x = 200, label.y = 0.052)+
  stat_cor(aes(label = after_stat(rr.label)), label.x = 200, label.y = 0.054)+
  scale_color_manual(name = "Population", values = mypalette) + 
  scale_fill_manual(name = "Population", values = mypalette) + 
  labs(title = "IBD in Arctic Cod",
       x = "Marine Distance (km)",
       y = expression("F"[ST]*"/(1-F"[ST]*")")) +
  theme_bw() + 
  theme(panel.border = element_blank(), axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 12),legend.title = element_text(size = 14))
totalIBD2
ggsave("./figures/ibd/boreogadus_allPops_ibd_r_r2.jpeg",
       totalIBD2, width = 8, height = 6, dpi = 300)

##### REMOVE ICELAND ##################
df.ibd2 <- df.ibd %>% filter(POP1!='Iceland',POP2!='Iceland')

westernIBD <- ggplot(df.ibd2, aes(x = geoDist, y = geneDist)) +
  geom_smooth(method = 'lm', se = F, col = "gray30") + 
  geom_point(size = 4.5, pch = 21, stroke = 2, aes(fill = POP1, color = POP2)) +
  #stat_cor(aes(label = after_stat(r.label)), label.x = 200, label.y = 0.052)+
  #stat_cor(aes(label = after_stat(rr.label)), label.x = 200, label.y = 0.054)+
  scale_color_manual(name = "Population", values = mypalette) + 
  scale_fill_manual(name = "Population", values = mypalette) + 
  labs(title = "IBD in Arctic Cod - Alaska and Canada",
       x = "Marine Distance (km)",
       y = expression("F"[ST]*"/(1-F"[ST]*")")) +
  theme_bw() + 
  theme(panel.border = element_blank(), axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 12),legend.title = element_text(size = 14))
westernIBD
ggsave("./figures/ibd/boreogadus_westernPops_ibd_colorByPop.jpeg",
       westernIBD, width = 8, height = 6, dpi = 300)

# all points one color
westernIBD2 <- ggplot(df.ibd2, aes(x = geoDist, y = geneDist)) +
  geom_smooth(method = 'lm', se = F, col = "gray60") + 
  geom_point(size = 4, pch = 21, fill = "darkslategray4", color = "black") +
  #stat_cor(aes(label = after_stat(r.label)), label.x = 200, label.y = 0.052)+
  #stat_cor(aes(label = after_stat(rr.label)), label.x = 200, label.y = 0.054)+
  scale_color_manual(name = "Population", values = mypalette) + 
  scale_fill_manual(name = "Population", values = mypalette) + 
  labs(title = "IBD in Arctic Cod - Alaska and Canada",
       x = "Marine Distance (km)",
       y = expression("F"[ST]*"/(1-F"[ST]*")")) +
  theme_bw() + 
  theme(panel.border = element_blank(), axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 12),legend.title = element_text(size = 14))
westernIBD2
ggsave("./figures/ibd/boreogadus_westernPops_ibd.jpeg",
       westernIBD2, width = 8, height = 6, dpi = 300)
