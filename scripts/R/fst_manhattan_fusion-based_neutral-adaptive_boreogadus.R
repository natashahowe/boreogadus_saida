# Boreogadus saida 
# z-score for Fst

### Install Packages & Load Functions ######################################

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "tools", "here", 
                     "NCmisc", "stringr") # these two are for z-scores

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

## Filepaths ################################################################
BAMFILE <- "./data/R/boreogadus_filtered_downsampled_bams.txt"
METADATAFILE <- "./data/R/boreogadus_metadata.csv"
COLORFILE <- paste0(here(),"./data/R/color_metadata_allpops_flip.txt")

### Metadata ###############################################################
chrom_df <- read.table("./data/R/Arctic_cod_genome_chromosomes.txt", header = T) # Read in tab-delimited table that has two columns chrom name from angsd and a simplified name (e.g., chr_1)

# simplified metadata
meta_df <- read.csv(METADATAFILE, header = T) %>%
  mutate(Population = ifelse(Population=='Southhampton Island', 'Southampton Island', Population),
         Pop = case_when(Population == 'SE Baffin' ~ 'SEBaffin',
                         TRUE ~ gsub(' .*$','',Population))) %>%
  filter(Pop != 'Unknown')

# create color palette
color_df <- read.delim2(COLORFILE, header = T, row.names = NULL, sep = "\t") %>%
  mutate(Pop = gsub('hh','h',Pop)) %>% # change Southhampton spelling
  left_join(meta_df %>% distinct(Population, Pop))

mypalette <- as.vector(color_df$Color) # turn colors into vector
  names(mypalette) <- color_df$Pop # attach pop name to palette color

POPLIST <- unique(meta_df$Pop)
  POPLIST
  POPLIST <- POPLIST[-grep("Labrador",POPLIST)]

fused_chroms <- c(1:5)
unfused_chroms <- c(6:18)

rm(COLORFILE, BAMFILE, METADATAFILE, meta_df)
###  Start for-loop by Pop ##################################################

i=9

# for(i in 5:length(POPLIST)){  
  
  FOCALPOP <- POPLIST[i]
  FOCALPOP
  
  # if focal pop = Southampton, change spelling to call in file with 'hh'
  FOCALPOP <- sub('Southampton','Southhampton',FOCALPOP) 
    
  WILDCARD <- paste0("*",FOCALPOP,"*.txt") # all text that include the name of the focal pop
  FILENAMES <- Sys.glob(file.path(here::here(),"results","fst",DATE,WILDCARD)) 
  basename(FILENAMES)
  if(length(grep('subset', FILENAMES))>0){
    FILENAMES <- FILENAMES[-grep("_subset",FILENAMES)] # not going to plot the subset pops in this
  }
  if(length(grep('Labrador', FILENAMES))>0){
    FILENAMES <- FILENAMES[-grep("Labrador",FILENAMES)]# error with labrador so not going to plot this round
  }
  basename(FILENAMES)
  #FILENAMES <- FILENAMES[1:3] #for testing
  
  #### Read in pairwise comparison data files ################################
  fst_df <- as.list(FILENAMES) %>%
    set_names(nm = FILENAMES) %>%
    map_dfr(
      ~ read_delim(.x, skip = 1, col_types = cols(), 
                   col_names = c("region", "chrName", "midPos", "Nsites", "Fst"), 
                   delim = "\t"), .id = "comparison"
    )
  
  rm(FILENAMES, WILDCARD)
  
  # edit dataframe
  fst_df <- inner_join(fst_df, chrom_df, by = "chrName") %>%
    select(chrName, chr, midPos, Fst, comparison) %>%
    mutate(midPos = as.numeric(midPos), Fst = as.numeric(Fst),
           temp = basename(comparison),
           POP_PAIR = str_extract(temp, "(?<=boreogadus_)[^_]+"),
           POP2 = str_replace(POP_PAIR, FOCALPOP, ""),
           POP2 = str_replace_all(POP2, "-", "")) %>%
    select(-c(comparison,temp,POP_PAIR))
  
  head(fst_df)
  
  # have to change POP2 name here as well
  fst_df <- fst_df %>%
    mutate(POP2 = ifelse(POP2=='Southhampton','Southampton',POP2),
           Fst = ifelse(Fst < 0, 0, Fst)) # change negative fst values to zero
  
  # Make the POP2 and chr column in the data frame into a factor with a specific order
  fst_df$POP2 <- factor(fst_df$POP2, levels = color_df$Pop)
  fst_df$chr <- factor(fst_df$chr, levels = chrom_df$chr)
  
  unique(fst_df$POP2)
  
  ### Fused: Adaptive Z-score calculations #######################################
  
  FusedPop <- fst_df %>%
    filter(chr %in% fused_chroms) %>%
    group_by(POP2) %>%                  # group data by Pop
    mutate(zFst = scale(Fst)[,1]) %>%   # z-statistic for Fst 
    ungroup()
  
  # Split by Pop and calculate n, p-val and z-score cutoff
  fused_zscores <- FusedPop %>%
    group_by(POP2) %>%                  # group data by Pop
    summarise(n = n(), avgFst = mean(Fst), sdFst = sd(Fst),
              min = min(zFst), max = max(zFst)) %>%              # calculate number of SNPs for each pop
    ungroup() %>%
    mutate(bonf_pval = alpha/n,
           Fst_1SD = p.to.Z(bonf_pval)*sdFst + avgFst)  # calculate the z-score cutoff value by Pop
  
  fused_zscores
  
  FusedPop <- FusedPop %>%
    left_join((fused_zscores %>% select(POP2,Fst_1SD)))
  
  # create dataframe to plot FST line (for each fused chr)
  FusedLine <- FusedPop %>% 
    distinct(POP2, chr, Fst_1SD)
  
  FusedLine$POP2 <- factor(FusedLine$POP2, levels = color_df$Pop)
  
  ### Unfused: Neutral Z-score calculations ######################################
  
  UnfusedPop <- fst_df %>%
    filter(chr %in% unfused_chroms) %>%
    group_by(POP2) %>%                  # group data by Pop
    mutate(zFst = scale(Fst)[,1]) %>%   # z-statistic for Fst
    ungroup()
  
  # Split by Pop and calculate n, p-val and z-score cutoff
  unfused_zscores <- UnfusedPop %>%
    group_by(POP2) %>%                  # group data by Pop
    summarise(n = n(), avgFst = mean(Fst), sdFst = sd(Fst),
              min = min(zFst), max = max(zFst)) %>%     # calculate number of SNPs for each pop
    ungroup() %>%
    mutate(bonf_pval = alpha/n,
           Fst_1SD = p.to.Z(bonf_pval)*sdFst + avgFst)  # calculate the z-score cutoff value by Pop
  
  unfused_zscores
  
  UnfusedPop <- UnfusedPop %>%
    left_join((unfused_zscores %>% select(POP2,Fst_1SD)))
  
  # create dataframe to plot FST line (for each unfused chr)
  UnfusedLine <- UnfusedPop %>% 
    distinct(POP2, chr, Fst_1SD)
  
  UnfusedLine$POP2 <- factor(UnfusedLine$POP2, levels = color_df$Pop)
  
  ##### Fst Plot ##############################################################
  
  # elongate the title to the Name of region rather than shortened pop name
  if(FOCALPOP=='Southhampton'){
    TITLE <- 'Southampton Island'
  }else{
    TITLE <- color_df$Population[which(color_df$Pop == FOCALPOP)]
  }
  TITLE
  
  # remove Fst=0 for plotting
  fst_df <- filter(fst_df,Fst > 0)
  
  # Set the ggplot theme
  theme_set(
    theme( 
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 70, size = 11, vjust = 0.5, color = "black"),
      axis.title.x = element_text(size = 16, color = "black"),
      axis.text.y = element_text(angle = 0, size = 11, color = "black"),
      axis.title.y = element_text(size = 20, color = "black", angle = 90),
      title = element_text(size = 15, color = "black"),
      panel.background = element_rect(fill = "white"), 
      panel.spacing = unit(0,"lines"),
      strip.text.x = element_text(angle = 0, color = "black", size = 14),
      strip.text.y = element_text(angle = 0, color = "black", size = 13),
      strip.background = element_rect(fill = "gray95")
    )
  )
  
  # Save the plot to a jpg file
  jpeg(paste0('./figures/fst/fusion/',FOCALPOP,'_fst_fusion_wholegenome_alpha',alpha,'_',format(Sys.Date(),format="%Y%m%d"),'.jpg'),
       width = 22, height = 10, res = 100, units = "in")

  print(ggplot() +
    geom_point(data = fst_df, aes(x = midPos/1e6, y = Fst, color = POP2), 
               alpha = 0.9, size = 0.9) +
    geom_hline(data = UnfusedLine, aes(yintercept = Fst_1SD), color = 'black', size = 1, linetype = "dashed") +
    geom_hline(data = FusedLine, aes(yintercept = Fst_1SD), color = 'black', size = 1, linetype = "solid") +
    scale_color_manual(values = mypalette) +
    facet_grid(POP2 ~ chr, scales = "free_x", space = "free_x",
               labeller = labeller(facet_category = label_wrap_gen(width = 13))) +
    ylab(expression(italic(F)[ST])) +
    xlab("Position (Mb)") +
    ggtitle(TITLE) +
    scale_y_continuous(breaks = seq(0.2, 1, by = 0.2)) +
    scale_x_continuous(breaks = seq(10,40, by = 10))
  )
  
  dev.off()
  
  rm(fst_df, FOCALPOP, TITLE,
     FusedLine, FusedPop, fused_zscores,  
     UnfusedLine, UnfusedPop, unfused_zscores)
#}  

# END OF GENOME PLOT LOOP


##############################################################################
# 
# # Save the plot to a jpg file
# jpeg(paste0('./figures/fst/zfst/',DATE,'/pval/',FOCALPOP,'_logp_wholegenome_alpha0.05-0.001_',format(Sys.Date(),format="%Y%m%d"),'.jpg'),
#      width = 22, height = 10, res = 100, units = "in")
# 
# ggplot() +
#   geom_point(data = zfst_df, aes(x = midPos/1e6, y = -log(pFst), color = POP2), 
#              alpha = 0.9, size = 0.9) +
#   geom_hline(data = zFst_stats, aes(yintercept = -log(bonf_pval*0.02)), color = 'black', linetype = "dashed") +
#   geom_hline(data = zFst_stats, aes(yintercept = -log(bonf_pval)), color = 'black', linetype = "dotted") +
#   scale_color_manual(values = mypalette) +
#   facet_grid(POP2 ~ chr, scales = "free", space = "free_x",
#              labeller = labeller(facet_category = label_wrap_gen(width = 13))) +
#   ylab('-log(p-value)') +
#   xlab("Position (Mb)") +
#   ggtitle(TITLE) +
#   scale_x_continuous(breaks = seq(10,40, by = 10))
# 
# dev.off()
# 
