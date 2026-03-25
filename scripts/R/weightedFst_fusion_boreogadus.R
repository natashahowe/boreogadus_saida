# This code calculated weighted FST separate for fused and unfused chroms
#
# If this has already been calculated and written to "./results/fst/boreogadus_fusion_weightedFSTs.csv",
# then can skip to section titled: Plotting Heatmap


packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "tidyverse",
                     "stringr", "data.table", "plyr","gtools", "tools","reshape2", 
                     "here", "magrittr", "patchwork", "cowplot")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

DATE="20241222" 

### Read in Metadata ####################################################

chrom_df <- read.table("./data/R/Arctic_cod_genome_chromosomes.txt", header = T) # Read in tab-delimited table that has two columns chrom name from angsd and a simplified name (e.g., chr_1)

meta_df <- read.delim("./data/R/region_metadata.txt", header = T) %>% # Read in metadata file with population, sampling location, and region
  mutate(Pop = sub('hh','h',Pop)) %>%
  filter(Pop != 'Labrador') # remove labrador
meta_df$Pop <- factor(meta_df$Pop, levels = unique(meta_df$Pop)) # Specify the order of some factors for plotting later

POPLIST <- unique(meta_df$Pop)
POPLIST

fused_chroms <- c(1:5)
unfused_chroms <- c(6:18)

### Import files for combining ##################################################
WEIGHTEDNAMES <- Sys.glob(file.path(here::here(),"results","fst",DATE,"print","*")) 
  #WEIGHTEDNAMES <- WEIGHTEDNAMES[-grep("_subset",WEIGHTEDNAMES)] # not going to plot the subset pops in this
  #WEIGHTEDNAMES <- WEIGHTEDNAMES[-grep("Labrador",WEIGHTEDNAMES)]# too few individuals in labrador so not going to plot this round
basename(WEIGHTEDNAMES)

# Read in pairwise comparison data files
idx_df <- as.list(WEIGHTEDNAMES) %>%
  set_names(nm = WEIGHTEDNAMES) %>%
  map_dfr(
    ~ read_delim(.x, col_types = cols(), 
                 col_names = c("chrName", "midPos", "A", "B"), delim = "\t"), .id = "filename"
  )
head(idx_df)

##### Combine each IDX Fst Comparison #######################
# bind rows together
idx_df <- idx_df %>%
  mutate(comparison = str_split(filename, "_", simplify = TRUE)[,3], 
         A = as.numeric(A), B = as.numeric(B),
         A = ifelse(A < 0, 0, A)) %>%
  left_join(chrom_df) %>%
  select(-filename)
  
unique(idx_df$comparison)

unfused_weightedFst <- idx_df %>%
  filter(chr %in% unfused_chroms) %>%
  group_by(comparison) %>%
  dplyr::summarise(unfusedFst = sum(A)/sum(B))
  
fused_weightedFst <- idx_df %>%
  filter(chr %in% fused_chroms) %>%
  group_by(comparison) %>%
  dplyr::summarise(fusedFst = sum(A)/sum(B))

rm(idx_df)

fusion_weightedFst <- inner_join(fused_weightedFst, unfused_weightedFst) %>%
  mutate(POP1 = sub("-.*","",comparison),
         POP2 = str_split(comparison, "-", simplify = TRUE)[,2])

write.csv(fusion_weightedFst, "./results/fst/boreogadus_fusion_weightedFst_dataframe.csv",
          quote = F, row.names = F)

#### Convert Dataframe to Matrix ###############

fusion_weightedFst <- read.csv("./results/fst/boreogadus_fusion_weightedFst_dataframe.csv",
                               row.names = NULL)

# have to reorder Pop1 and Pop2 data to reflect location order
fusion_weightedFst <- fusion_weightedFst %>%
  mutate(Pop1 = case_when(match(POP1, levels(POPLIST)) < match(POP2, levels(POPLIST)) ~ POP1,
                          TRUE ~ POP2),
         Pop2 = case_when(Pop1 == POP1 ~ POP2, 
                          TRUE ~ POP1)) %>%
  select(-c(comparison, POP1, POP2))

# create row and col names and order appropriately
rows <- levels(POPLIST)[(levels(POPLIST) %in% unique(fusion_weightedFst$Pop1))]
cols <- levels(POPLIST)[(levels(POPLIST) %in% unique(fusion_weightedFst$Pop2))]

# Make sure rows and cols are the same and in the same order
all_names <- unique(c(rows, cols))

# Initialize a square zero matrix
mat <- matrix(0, nrow = length(all_names), ncol = length(all_names),
              dimnames = list(all_names, all_names))

# Now fill the matrix
for (i in 1:nrow(fusion_weightedFst)) {
  row_name <- fusion_weightedFst$Pop1[i]
  col_name <- fusion_weightedFst$Pop2[i]
  fused_Fst <- fusion_weightedFst$fusedFst[i]
  unfused_Fst <- fusion_weightedFst$unfusedFst[i]
  
  row_idx <- match(row_name, all_names)
  col_idx <- match(col_name, all_names)
  
  if (row_idx > col_idx) {
    mat[row_idx, col_idx] <- fused_Fst  # Lower triangle
    mat[col_idx, row_idx] <- unfused_Fst  # Upper triangle
  } else if (row_idx < col_idx) {
    mat[row_idx, col_idx] <- unfused_Fst  # Upper triangle
    mat[col_idx, row_idx] <- fused_Fst  # Lower triangle
  }
  # Diagonal remains 0
}

mat

write.csv(mat, "./results/fst/boreogadus_fusion_weightedFSTs.csv",
            row.names = T, quote = F)

### Plotting Heatmap ###########################################################
mat <- read.csv("./results/fst/boreogadus_fusion_weightedFSTs.csv")

rownames(mat) <- mat[,1]
mat <- as.matrix(mat[,2:ncol(mat)])

df_long <- as.data.frame(as.table(mat))
colnames(df_long) <- c("Pop1", "Pop2", "Fst")
df_long$Fst[df_long$Fst == 0] <- NA

df_long <- df_long %>%
  mutate(Pop1 = gsub('hh','h',Pop1),
         Pop2 = gsub('hh','h',Pop2))

df_long <- df_long %>%
  mutate(Pop1 = factor(Pop1, levels = POPLIST),
         Pop2 = factor(Pop2, levels = POPLIST))

heatmap <- ggplot(df_long, aes(Pop1, Pop2, fill = Fst)) +
  geom_tile(color = "white") +
  geom_text(aes(label = format(round(Fst, 3), nsmall = 2)), size = 3.5, na.rm = T) +  
  scale_fill_distiller(palette = "Reds", direction = 1,
                       na.value = "white", name = expression(italic(F)[ST])) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11)) +
  coord_fixed() + labs(x = "", y = "")
heatmap

# ggsave("./figures/fst/fusion/boreogadus_fusion_weightedFst_heatmap.jpeg",
#        heatmap, width = 8, height = 7)

# different color scale
heatmap2 <- heatmap +  
  scale_fill_distiller(palette = "Spectral", na.value = "white", name = expression(italic(F)[ST]))
heatmap2

ggsave("./figures/fst/fusion/boreogadus_fusion_weightedFst_heatmap_spectral.jpeg",
       heatmap2, width = 8, height = 7)
