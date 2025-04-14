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

##  Start For Loop by Pop ##################################################
# identify files for combining
  Adaptive_Pop <- list()

i=7

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
  
  head(fst_df)
  unique(fst_df$POP2)
  
  # Make the POP2 and chr column in the data frame into a factor with a specific order
  fst_df$POP2 <- factor(fst_df$POP2, levels = unique(meta_df$Pop))
  fst_df$chr <- factor(fst_df$chr, levels = chrom_df$chr)
  
  # get rid of negative fst values
  fst_df$Fst[fst_df$Fst < 0] <- 0
  
  #### Z-score calculations for Outliers ############

  # Split by Pop and calculate n, p-val and z-score cutoff
  zscores <- fst_df %>%
    group_by(POP2) %>%                  # Split up following calcs by Pop
    summarise(n = n()) %>%              # calculate number of SNPs for each pop
    ungroup() %>%                       # now it is split by pop with number
    mutate(bonf_pval = alpha/n)  # calculate the z-score cutoff value by Pop
  
  zfst_df <- fst_df %>%
    group_by(POP2) %>%                  # group data by Pop
    mutate(zFst = scale(Fst)[,1],       # z-statistic for Fst
           pFst = 2*pnorm(-abs(zFst), mean=0, sd=2)) %>% # p-val with 2SD
    ungroup()
  
  # View min and max of z-scores
  zFst_stats <- zfst_df %>%
    group_by(POP2) %>%                  # group data by Pop
    summarise(avgFst = mean(Fst), sdFst = sd(Fst),
              min = min(zFst), max = max(zFst)) %>%
    ungroup() %>%
    left_join(zscores, by = "POP2") %>%                # combine with cutoff data
    mutate(Fst_2SD = (abs(qnorm(bonf_pval/ 2))*2)*sdFst + avgFst,
           Fst_1SD = p.to.Z(bonf_pval)*sdFst + avgFst) # reverse calculate FST from z-score = (FST - meanFST)/sdFST
  zFst_stats
  
  rm(fst_df, zscores)
  
  ##### Fst Plot ##############################################################
  
  # elongate the title to the Name of region rather than shortened pop name
  TITLE <- meta_df$Name[which(meta_df$Pop == FOCALPOP)]
  
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
      strip.background = element_rect(fill = "gray90")
    )
  )
  
  # Save the plot to a jpg file
  jpeg(paste0('./figures/fst/zfst/',DATE,'/2SD/',FOCALPOP,'_fst_wholegenome_alpha',alpha,'_',format(Sys.Date(),format="%Y%m%d"),'.jpg'),
       width = 22, height = 10, res = 100, units = "in")

  ggplot() +
    geom_point(data = zfst_df, aes(x = midPos/1e6, y = Fst, color = POP2), 
               alpha = 0.9, size = 0.9) +
    geom_hline(data = zFst_stats, aes(yintercept = Fst_1SD), color = 'black', linetype = "dashed") +
    geom_hline(data = zFst_stats, aes(yintercept = Fst_2SD), color = 'black', linetype = "dotted") +
    scale_color_manual(values = mypalette) +
    facet_grid(POP2 ~ chr, scales = "free_x", space = "free_x",
               labeller = labeller(facet_category = label_wrap_gen(width = 13))) +
    ylab(expression(italic(F[ST]))) +
    xlab("Position (Mb)") +
    ggtitle(TITLE) +
    scale_y_continuous(breaks = seq(0.2, 1, by = 0.2)) +
    scale_x_continuous(breaks = seq(10,40, by = 10))
  
  dev.off()
  
  rm(zfst_df, zFst_stats, FOCALPOP)
}  

# END OF GENOME PLOT LOOP


##############################################################################

# Save the plot to a jpg file
jpeg(paste0('./figures/fst/zfst/',DATE,'/pval/',FOCALPOP,'_logp_wholegenome_alpha0.05-0.001_',format(Sys.Date(),format="%Y%m%d"),'.jpg'),
     width = 22, height = 10, res = 100, units = "in")

ggplot() +
  geom_point(data = zfst_df, aes(x = midPos/1e6, y = -log(pFst), color = POP2), 
             alpha = 0.9, size = 0.9) +
  geom_hline(data = zFst_stats, aes(yintercept = -log(bonf_pval*0.02)), color = 'black', linetype = "dashed") +
  geom_hline(data = zFst_stats, aes(yintercept = -log(bonf_pval)), color = 'black', linetype = "dotted") +
  scale_color_manual(values = mypalette) +
  facet_grid(POP2 ~ chr, scales = "free", space = "free_x",
             labeller = labeller(facet_category = label_wrap_gen(width = 13))) +
  ylab('-log(p-value)') +
  xlab("Position (Mb)") +
  ggtitle(TITLE) +
  scale_x_continuous(breaks = seq(10,40, by = 10))

dev.off()

