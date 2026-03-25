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
color_df <- read.delim("./data/R/color_metadata_allpops_flip.txt", header = T) %>%
  mutate(Pop = gsub('hh','h',Pop)) # correct spelling

POPLIST = color_df$Pop
POPLIST

# create mypalette
mypalette <- as.vector(color_df$Color) # turn colors into vector
  names(mypalette) <- POPLIST # attach pop name to palette color
  mypalette

#################################################################################
# files for combining
  WEIGHTEDNAMES <- Sys.glob(file.path(here::here(),"results","fst",DATE,"global","*20241222*")) 
    WEIGHTEDNAMES <- WEIGHTEDNAMES[-grep("_subset",WEIGHTEDNAMES)] # not going to plot the subset pops in this
    WEIGHTEDNAMES <- WEIGHTEDNAMES[-grep("Labrador",WEIGHTEDNAMES)]# too few individuals in labrador so not going to plot this round
  basename(WEIGHTEDNAMES)
  weighted_list <- as.list(WEIGHTEDNAMES)
  
  ################################################################################################
  # Read in pairwise comparison data files
  weighted_fst <- weighted_list %>%
    set_names(nm = WEIGHTEDNAMES) %>%
    map_dfr(
      ~ read_delim(.x, col_types = cols(), 
                   col_names = c("unweighted", "weighted_fst"), delim = "\t"), .id = "filename"
    )
  head(weighted_fst)
  
  # Strip the original file name such that only the population that FOCALPOP is 
  # being compared to is retained (for plot labeling)
  weighted_fst <- weighted_fst %>%
    mutate(comparison = str_split(filename, "_", simplify = TRUE)[,3],
           POP1 = gsub("-.*","",comparison),
           POP2 = gsub(".*-","",comparison)) %>%
    select(c(POP1, POP2, weighted_fst)) %>%
    mutate(POP1 = gsub('hh','h',POP1),
           POP2 = gsub('hh','h',POP2))
  
  head(weighted_fst)
  unique(weighted_fst$POP2)

  # Make the POP2 column in the data frame into a factor with a specific order
  weighted_fst$POP1 <- factor(weighted_fst$POP1, levels = POPLIST)
  weighted_fst$POP2 <- factor(weighted_fst$POP2, levels = POPLIST)
  
  weighted_fst_rev <- weighted_fst %>%
    dplyr::rename(POP2 = POP1, POP1 = POP2) %>% select(POP1,everything())

  mat.fst <- acast(rbind(weighted_fst,weighted_fst_rev), POP1~POP2, value.var = "weighted_fst")  

## Specify some functions to retrieve upper part of matrix
# Get lower triangle of the correlation matrix

get_lower_tri <- function(Fstmat){
  Fstmat[upper.tri(Fstmat)] <- NA
  return(Fstmat)
}

## subset the matrix
lower_tri <- get_lower_tri(mat.fst)
View(lower_tri)

##Use the package reshape to melt the matrix into a df again:
final_df <- melt(lower_tri, value.name = "weighted_fst") %>%
  filter(!is.na(weighted_fst)) %>%
  mutate(weighted_fst = round(weighted_fst, digits = 4))

# order populations by the color_df dataframe order
final_df$POP1 <- factor(as.factor(final_df$Var1), levels = POPLIST)
final_df$POP2 <- factor(final_df$Var2, levels = POPLIST)

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
  ylab("") + xlab("") +
  labs(fill = expression(italic(F)[ST])) +
  coord_fixed()
heatmap_plot

# save plot to file
ggsave(paste0("./figures/fst/",DATE,"/weighted_Fst_heatmap2.jpeg"), 
       heatmap_plot, width= 7, height = 6)

#### ADD FUSION DIFFERENCE MATRIX ################################################

fusion_weightedFst <- read.csv("./results/fst/boreogadus_fusion_weightedFst_dataframe.csv",
                               row.names = NULL) %>%
  mutate(POP1 = gsub('hh','h',POP1), POP2 = gsub('hh','h',POP2))

# have to reorder Pop1 and Pop2 data to reflect location order
fusion_Fst <- fusion_weightedFst %>%
  mutate(diffFst = fusedFst - unfusedFst,
         propFst = fusedFst / unfusedFst, 
         Pop1 = case_when(match(POP1, POPLIST) < match(POP2, POPLIST) ~ POP1,
                          TRUE ~ POP2),
         Pop2 = case_when(Pop1 == POP1 ~ POP2, 
                          TRUE ~ POP1)) %>%
  select(-POP1, -POP2, -comparison) 

diff.global.fst <- fusion_Fst %>%
  left_join(weighted_fst, by = join_by("Pop1" == "POP1", "Pop2" == "POP2")) %>%
  left_join(weighted_fst, by = join_by("Pop2" == "POP1", "Pop1" == "POP2")) %>%
  mutate(globalFst = ifelse(is.na(weighted_fst.x), weighted_fst.y, weighted_fst.x)) %>%
  select(-starts_with("weight"))

# create row and col names and order appropriately
rows <- POPLIST[(POPLIST %in% unique(diff.global.fst$Pop1))]
cols <- POPLIST[(POPLIST %in% unique(diff.global.fst$Pop2))]

# Make sure rows and cols are the same and in the same order
all_names <- unique(c(rows, cols))
all_names

# Initialize a square zero matrix
diff.global.mat <- matrix(0, nrow = length(all_names), ncol = length(all_names),
              dimnames = list(all_names, all_names))

# Now fill the matrix
for (i in 1:nrow(diff.global.fst)) {
  row_name <- diff.global.fst$Pop1[i]
  col_name <- diff.global.fst$Pop2[i]
  weighted_Fst <- diff.global.fst$globalFst[i]
  diff_Fst <- diff.global.fst$diffFst[i]
  
  row_idx <- match(row_name, all_names)
  col_idx <- match(col_name, all_names)
  
  if (row_idx > col_idx) {
    diff.global.mat[row_idx, col_idx] <- diff_Fst  # Lower triangle
    diff.global.mat[col_idx, row_idx] <- weighted_Fst  # Upper triangle
  } else if (row_idx < col_idx) {
    diff.global.mat[row_idx, col_idx] <- weighted_Fst  # Upper triangle
    diff.global.mat[col_idx, row_idx] <- diff_Fst  # Lower triangle
  }
  # Diagonal remains 0
}

diff.global.mat

## plot 
plot_df <- as.data.frame(as.table(diff.global.mat))
colnames(plot_df) <- c("Pop1", "Pop2", "Fst")
plot_df$Fst[plot_df$Fst == 0] <- NA

diff.global.heatmap <- ggplot(plot_df, aes(Pop1, Pop2, fill = Fst)) +
  geom_tile(color = "white") +
  geom_text(aes(label = format(round(Fst, 3), nsmall = 2)), size = 3.5, na.rm = T) +  
  scale_fill_distiller(palette = "Spectral", 
                       na.value = "white", name = expression(italic(F)[ST])) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11)) +
  coord_fixed() + labs(x = "", y = "")
diff.global.heatmap

ggsave("./figures/fst/fusion/boreogadus_wholegenomeFst_and_fusionDiffFst_heatmap_spectral.jpeg",
       diff.global.heatmap, width = 8, height = 7)

########## PROPORTIONAL FST #####################################

# Initialize a square zero matrix
prop.global.mat <- matrix(0, nrow = length(all_names), ncol = length(all_names),
                          dimnames = list(all_names, all_names))

# Now fill the matrix
for(i in 1:nrow(diff.global.fst)) {
  row_name <- diff.global.fst$Pop1[i]
  col_name <- diff.global.fst$Pop2[i]
  weighted_Fst <- diff.global.fst$globalFst[i]
  prop_Fst <- diff.global.fst$propFst[i]
  
  row_idx <- match(row_name, all_names)
  col_idx <- match(col_name, all_names)
  
  if (row_idx > col_idx) {
    prop.global.mat[row_idx, col_idx] <- weighted_Fst  # Lower triangle
    prop.global.mat[col_idx, row_idx] <- prop_Fst  # Upper triangle
  } else if (row_idx < col_idx) {
    prop.global.mat[row_idx, col_idx] <- prop_Fst  # Upper triangle
    prop.global.mat[col_idx, row_idx] <- weighted_Fst  # Lower triangle
  }
  # Diagonal remains 0
}

prop.global.mat

## plot 
prop.plot_df <- as.data.frame(as.table(prop.global.mat))
colnames(prop.plot_df) <- c("Pop1", "Pop2", "Fst")
prop.plot_df$Fst[prop.plot_df$Fst == 0] <- NA

prop.global.heatmap <- ggplot(prop.plot_df, aes(Pop1, Pop2, fill = Fst)) +
  geom_tile(color = "white") +
  geom_text(aes(label = format(round(Fst, 3), nsmall = 2)), size = 3.5, na.rm = T) +  
  scale_fill_distiller(palette = "Spectral", 
                       na.value = "white", name = expression(italic(F)[ST])) +
theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11)) +
  coord_fixed() + labs(x = "", y = "")
prop.global.heatmap

ggsave("./figures/fst/fusion/boreogadus_globalFst_and_fusionPropFst_heatmap_spectral.jpeg",
       prop.global.heatmap, width = 8, height = 7)
