# Boreogadus saida 
# z-score for Fst

######################################
### INSTALL PACKAGES & LOAD FUNCTIONS

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "tools", "here", 
                     "NCmisc", "scales") # these two are for z-scores

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

# CHOOSE OPTIONS
# for folder directory to angsd output and plot output
# for this code to work, will also need to match date in file names
# for example: Broughton.Frobisher_20241222.fst.SNP.txt
DATE="20241222" 

#################################################################################
# Read in tab-delimited table that has two columns chrom name from angsd and a simplified name (e.g., chr_1)
chrom_df <- read.table("./data/R/Arctic_cod_genome_chromosomes.txt", header = T)

# Read in metadata file with population, sampling location, and region
meta_df <- read.delim("./data/R/region_metadata.txt", header = T) 

# read in color file
color_df <- read.delim("./data/R/color_metadata_allpops_flip.txt", header = T) 

# Specify the order of some factors for plotting later
meta_df$Pop <- factor(meta_df$Pop, levels = meta_df$Pop)
meta_df$Region <- factor(meta_df$Region, levels = unique(meta_df$Region))

POPLIST <- unique(meta_df$Pop)
POPLIST

# create mypalette
mypalette <- as.vector(color_df$Color) # turn colors into vector
names(mypalette) <- color_df$Pop # attach pop name to palette color
mypalette


#################################################################################
# identify files for combining
i=3

for(i in 1:length(meta_df$Pop)){
  
  FOCALPOP <- as.character(meta_df$Pop[i])
  FOCALPOP
  
  WILDCARD <- paste0("*",FOCALPOP,"*.txt") # all text that include the name of the focal pop
  FILENAMES <- Sys.glob(file.path(here(),"results","fst",DATE,WILDCARD)) 
  basename(FILENAMES)
  if(length(grep('subset', FILENAMES))>0){
    FILENAMES <- FILENAMES[-grep("_subset",FILENAMES)] # not going to plot the subset pops in this
  }
  if(length(grep('Labrador', FILENAMES))>0){
    FILENAMES <- FILENAMES[-grep("Labrador",FILENAMES)]# error with labrador so not going to plot this round
  }
  basename(FILENAMES)
  file_list <- as.list(FILENAMES)
  
  ################################################################################################
  # Read in pairwise comparison data files
  fst_df <- file_list %>%
    set_names(nm = FILENAMES) %>%
    map_dfr(
      ~ read_delim(.x, skip = 1, col_types = cols(), 
                   col_names = c("region", "chrName", "midPos", "Nsites", "Fst"), 
                   delim = "\t"), .id = "comparison"
    )
  
  # Edit data frame
  fst_df <- inner_join(fst_df, chrom_df, by = "chrName") %>%
    select(chrName, chr, midPos, Fst, comparison) %>%
    mutate(midPos = as.numeric(midPos)/1e6,
           Fst = as.numeric(Fst))
  
  # Strip the original file name such that only the population that FOCALPOP is 
  # being compared to is retained (for plot labeling)
  fst_df <- fst_df %>%
    mutate(temp = basename(tools::file_path_sans_ext(comparison)),
           temp = sub("boreogadus_","",temp),
           temp = sub("_subset","",temp),
           temp = sub("_downsampled","",temp),
           temp = gsub(paste0("_",DATE,".fst.SNP"), "", temp),
           temp = gsub(FOCALPOP, "", temp),
           POP2 = gsub("-", "", temp)) %>%
    select(-c(comparison, temp))
  
  # for z-score testing
  fst_df <- fst_df %>% filter(POP2 == "Broughton")
  
  head(fst_df)
  unique(fst_df$POP2)
  
  # Make the POP2 and chr column in the data frame into a factor with a specific order
  fst_df$POP2 <- factor(fst_df$POP2, levels = unique(meta_df$Pop))
  fst_df$chr <- factor(fst_df$chr, levels = chrom_df$chr)
  
  # get rid of negative fst values
  fst_df$Fst[fst_df$Fst < 0] <- 0
  
  #### Z-score calculations for Outliers ############
  library(NCmisc)
  
  alpha = 0.001
  
  # Split by Pop and calculate n, p-val and z-score cutoff
  zscores <- fst_df %>%
    group_by(POP2) %>%            # Split up following calcs by Pop
    summarise(n = n()) %>%        # calculate number of SNPs for each pop
    ungroup() %>%                 # now it is split by pop with number
    mutate(bonf_pval = alpha/n,   # calculate bonferonni corrected p-val
           zscore = p.to.Z(bonf_pval))  # calculate the z-score cutoff
  
  zfst_df <- fst_df %>%
    group_by(POP2) %>%
    mutate(zFst = scale(Fst)[,1]) %>%
    ungroup()
  
  write.table(zfst_df, "./results/Pond-Broughton_zFst.txt",
              row.names = F, col.names = T, quote = F)
  
  # View min and max of z-scores
  zFst_stats <- zfst_df %>%
    group_by(POP2) %>%
    summarise(avgFst = mean(Fst), sdFst = sd(Fst), 
              min = min(zFst), max = max(zFst)) %>%
    ungroup() %>%
    left_join(zscores, by = "POP2") %>%
    mutate(Fst_cutoff = zscore*sdFst + avgFst) # reverse calculating z-score = (FST - meanFST)/sdFST
  zFst_stats