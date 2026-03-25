# Find Adaptive & Neutral Loci for B saida

# Boreogadus saida 
# z-score for Fst

### Install Packages & Load Functions ######################################

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "tools", "here", 
                     "NCmisc", "scales") # these two are for z-scores

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}
rm(packages_needed, i)

##### CHOOSE OPTIONS ##########################################################
# for folder directory to angsd output and plot output
# for this code to work, will also need to match date in file names
# for example: Broughton.Frobisher_20241222.fst.SNP.txt

DATE="20241222" 
alpha = 0.05 # statistical alpha. Tested 0.05 and 0.001

### Read in Metadata ####################################################

chrom_df <- read.table("./data/R/Arctic_cod_genome_chromosomes.txt", header = T) # Read in tab-delimited table that has two columns chrom name from angsd and a simplified name (e.g., chr_1)

meta_df <- read.delim("./data/R/region_metadata.txt", header = T) %>% # Read in metadata file with population, sampling location, and region
  filter(Pop != 'Labrador') # remove labrador
  meta_df$Pop <- factor(meta_df$Pop, levels = meta_df$Pop) # Specify the order of some factors for plotting later
  meta_df$Region <- factor(meta_df$Region, levels = unique(meta_df$Region)) # Specify the order of some factors for plotting later

POPLIST <- unique(meta_df$Pop)
  POPLIST

color_df <- read.delim("./data/R/color_metadata_allpops_flip.txt", header = T) # read in color file
  mypalette <- as.vector(color_df$Color) # turn colors into vector
  names(mypalette) <- color_df$Pop # attach pop name to palette color
  mypalette

fused_chroms <- c(1:5)
unfused_chroms <- c(6:18)

##  Start For Loop by Pop ##################################################
# identify files for combining

FusedPop <- list(); UnfusedPop <- list()
AdaptivePop <- list(); NeutralPop <- list(); UnfusedAdaptivePop <- list()

#i=3

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
  
  ### TESTING SUBSET ####
  #FILENAMES <- FILENAMES[1:2]   ###### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  file_list <- as.list(FILENAMES)
  
  #### Read in pairwise comparison data files ################################
  fst_df <- file_list %>%
    set_names(nm = FILENAMES) %>%
    map_dfr(
      ~ read_delim(.x, skip = 1, col_types = cols(), 
                   col_names = c("region", "chrName", "midPos", "Nsites", "Fst"), 
                   delim = "\t"), .id = "comparison"
    )
  
  rm(file_list, FILENAMES, WILDCARD)
  
  # Edit data frame
  fst_df <- inner_join(fst_df, chrom_df, by = "chrName") %>%
    select(chrName, chr, midPos, Fst, comparison) %>%
    mutate(midPos = as.numeric(midPos),
           Fst = as.numeric(Fst),
           temp = basename(tools::file_path_sans_ext(comparison)),
           temp = sub("boreogadus_","",temp),
           temp = sub("_subset","",temp),
           temp = sub("_downsampled","",temp),
           temp = gsub(paste0("_",DATE,".fst.SNP"), "", temp),
           temp = gsub(FOCALPOP, "", temp),
           POP2 = gsub("-", "", temp)) %>%
    select(-c(comparison, temp))

fst_df$chr <- factor(fst_df$chr, levels = unique(fst_df$chr))

# get rid of negative fst values
fst_df$Fst[fst_df$Fst < 0] <- 0
  head(fst_df)
  unique(fst_df$POP2)

# split into fused and un-fused chromosomes
fused_df <- fst_df %>%
  filter(chr %in% fused_chroms)

unfused_df <- fst_df %>%
  filter(chr %in% unfused_chroms)

  unique(fused_df$chr)
  unique(unfused_df$chr)
  
  rm(fst_df)

### Fused: Adaptive Z-score calculations #######################################

FusedPop[[i]] <- fused_df %>%
  group_by(POP2) %>%                  # group data by Pop
  mutate(zFst = scale(Fst)[,1]) %>%   # z-statistic for Fst
  ungroup()

# Split by Pop and calculate n, p-val and z-score cutoff
fused_zscores <- FusedPop[[i]] %>%
  group_by(POP2) %>%                  # group data by Pop
  summarise(n = n(), avgFst = mean(Fst), sdFst = sd(Fst),
            min = min(zFst), max = max(zFst)) %>%              # calculate number of SNPs for each pop
  ungroup() %>%
  mutate(bonf_pval = alpha/n,
         Fst_1SD = p.to.Z(bonf_pval)*sdFst + avgFst)  # calculate the z-score cutoff value by Pop

fused_zscores

#### Fused Adaptive Loci #######################

Adaptive <- list()

for(j in 1:length(fused_zscores$POP2)){
  Adaptive[[j]] <- FusedPop[[i]] %>%
    filter(POP2 == fused_zscores$POP2[j],
           Fst > fused_zscores$Fst_1SD[j])
}

AdaptivePop[[i]] <- do.call(rbind, Adaptive) %>%
  mutate(markerPos = paste0(chrName,'_',midPos)) %>%
  select(markerPos) %>%
  distinct()
nrow(AdaptivePop[[i]])
# 
### Unfused: Neutral Z-score calculations ######################################

UnfusedPop[[i]] <- unfused_df %>%
  group_by(POP2) %>%                  # group data by Pop
  mutate(zFst = scale(Fst)[,1]) %>%   # z-statistic for Fst
  ungroup()

# Split by Pop and calculate n, p-val and z-score cutoff
unfused_zscores <- UnfusedPop[[i]] %>%
  group_by(POP2) %>%                  # group data by Pop
  summarise(n = n(), avgFst = mean(Fst), sdFst = sd(Fst),
            min = min(zFst), max = max(zFst)) %>%              # calculate number of SNPs for each pop
  ungroup() %>%
  mutate(bonf_pval = alpha/n,
         Fst_1SD = p.to.Z(bonf_pval)*sdFst + avgFst)  # calculate the z-score cutoff value by Pop

unfused_zscores

##### Neutral Loci ######################

# This is done by finding Unfused Adaptive < 1SD and removing them from the Neutral Dataset

UnfusedAdaptive <- list()

for(j in 1:length(unfused_zscores$POP2)){
  UnfusedAdaptive[[j]] <- UnfusedPop[[i]] %>%
    filter(POP2 == unfused_zscores$POP2[j],
           Fst > unfused_zscores$Fst_1SD[j])

}

UnfusedAdaptivePop[[i]] <- do.call(rbind, UnfusedAdaptive) %>%
  mutate(markerPos = paste0(chrName,'_',midPos)) %>% select(markerPos)

NeutralPop[[i]] <- unfused_df %>%
  mutate(markerPos = paste0(chrName,'_',midPos)) %>%
  filter(!(markerPos %in% UnfusedAdaptivePop[[i]]$markerPos)) %>% # remove if it is adaptive in any comparison
  select(markerPos) %>%
  distinct()

nrow(NeutralPop[[i]])

rm(FOCALPOP, Adaptive, fused_zscores, unfused_zscores, UnfusedAdaptive)

}  

### END OF GENOME PLOT LOOP


# # Combine Adaptive & Neutral SNPs across Pops
# adaptiveTotal_df <- do.call(rbind, AdaptivePop) %>%
#   arrange(markerPos) %>%
#   distinct(markerPos)
# 
# nrow(adaptiveTotal_df)  
# 
# write.table(adaptiveTotal_df, paste0("./results/fst/",DATE,"/NeutralAdaptiveMarkers/Fusion/boreogadus_fused_adaptive_markerPos_alpha",alpha,".txt"),
#             quote=F, row.names = F, col.names = F)

#### concatenate unfused adaptive markers

unfusedAdaptiveTotal_df <- do.call(rbind, UnfusedAdaptivePop) %>%
  arrange(markerPos) %>%
  distinct(markerPos)

tail(unfusedAdaptiveTotal_df)

write.table(unfusedAdaptiveTotal_df, paste0("./results/fst/",DATE,"/NeutralAdaptiveMarkers/Fusion/boreogadus_unfused_adaptive_markerPos_alpha",alpha,".txt"),
            quote=F, row.names = F, col.names = F)

# neutralTotal_df <- do.call(rbind, NeutralPop) %>%
#   arrange(markerPos) %>%
#   distinct(markerPos) %>%
#   filter(!(markerPos %in% unfusedAdaptiveTotal_df$markerPos)) # remove if it adaptive in any comparison
# 
# nrow(neutralTotal_df)
# 
# write.table(neutralTotal_df, paste0("./results/fst/",DATE,"/NeutralAdaptiveMarkers/Fusion/boreogadus_unfused_neutral_markerPos_alpha",alpha,".txt"),
#             quote=F, row.names = F, col.names = F)

