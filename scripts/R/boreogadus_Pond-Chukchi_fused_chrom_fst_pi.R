######################################
### INSTALL PACKAGES & LOAD FUNCTIONS

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "tools", 
                     "here", "patchwork")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

DATE="20241222" 

#################################################################################
  # FST file
  FSTFILE <- "./results/fst/20241222/boreogadus_Chukchi-Pond_20241222.fst.SNP.txt"

  # Read in tab-delimited table that has two columns chrom name from angsd and a simplified name (e.g., chr_1)
  chrom_df <- read.table("./data/R/Arctic_cod_genome_chromosomes.txt", header = T)

  # Read in metadata file with population, sampling location, and region
  meta_df <- read.delim("./data/R/region_metadata.txt", header = T)
  
  # read in color file
  color_df <- read.delim("./data/R/color_metadata_allpops_flip.txt", header = T)

  fusionPoint <- data.frame(chr = 1:5,
                            minPos = c(24,26,36,22,18.5),
                            maxPos = c(30,35,40,24,20))
  
  centromere <- data.frame(chr = 1:5,
                           minPos = c(29,30,5,19,13),
                           maxPos = c(30,31,6,20,14))
  
  ################################################################################################
  # Read in pairwise comparison data files
  fst_df <- read_delim(FSTFILE, skip = 1, 
                       col_types = cols(), col_names = c("region", "chrName", "midPos", "Nsites", "Fst"),
                       delim = "\t")

  # Edit data frame
  fst_df <- inner_join(fst_df, chrom_df, by = "chrName") %>%
    select(chr, midPos, Fst) %>%
    mutate(midPos = as.numeric(midPos)/1e6,
           Fst = as.numeric(Fst),
           Fst = ifelse(Fst < 0, 0, Fst))
  head(fst_df)
  
  # only keep 1-5 && add centromere locations
  fused_df <- fst_df %>%
    filter(chr <= 5) %>%
    rowwise() %>%
    mutate(fusionPoint = any(chr == fusionPoint$chr &
                               midPos >= fusionPoint$minPos &
                               midPos <= fusionPoint$maxPos),
           centromere = any(chr == centromere$chr &
                              midPos >= centromere$minPos &
                              midPos <= centromere$maxPos),
           CategorizeSNP = case_when(centromere == T ~ "centromere",
                                     fusionPoint == T ~ "fusion",
                                     T ~ "neither"))
  
  head(fused_df)
  unique(fused_df$chr)

  ##### Add in Pi #####################################################
  
  PIFILES <- Sys.glob(file.path(here(),"results","pi","pop_pi","b*")) 
  pi_dfs <- NULL
  for(i in 1:length(PIFILES)){
    pi_dfs[[i]] <- read_csv(PIFILES[i], col_types = cols())
  }
  
  head(pi_dfs[[1]])
  
  do.call(rbind,pi_dfs) %>%
    count(chrName)
  
  # bind all chromosomes
  pi_df <- do.call(rbind,pi_dfs) %>%
    mutate(midPos = as.numeric(midPos)/1e6)
  
  pi_df$Chukchi_pi[startsWith(pi_df$Chukchi_pi, "#")] <- NA
  pi_df$Pond_pi[startsWith(pi_df$Pond_pi, "#")] <- NA
  
  pi_df <- pi_df %>%
    mutate(Pond_pi = as.numeric(Pond_pi),
           Chukchi_pi = as.numeric(Chukchi_pi),
           delta_pi = Pond_pi - Chukchi_pi) %>%
    left_join(chrom_df) %>%
  # add fusionPoint locations
    rowwise() %>%
    mutate(fusionPoint = any(chr == fusionPoint$chr &
                               midPos >= fusionPoint$minPos &
                               midPos <= fusionPoint$maxPos),
           centromere = any(chr == centromere$chr &
                              midPos >= centromere$minPos &
                              midPos <= centromere$maxPos),
           CategorizeSNP = case_when(centromere == T ~ "centromere",
                                     fusionPoint == T ~ "fusion",
                                     T ~ "neither"))
  
  #### Part 2: Plot the data ###################################################
  
  # Set the ggplot theme
  theme_set(
    theme( 
      legend.position = "none",
      panel.grid.major = element_line(color = "gray95"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, color = "black"),
      axis.title.x = element_text(size = 16, color = "black"),
      axis.text.y = element_text(angle = 0, size = 12, color = "black"),
      axis.title.y = element_text(size = 20, color = "black", angle = 90),
      title = element_blank(),
      panel.background = element_rect(fill = "white"), 
      panel.spacing = unit(0,"lines"),
      strip.text.x = element_text(angle = 0, color = "black", size = 15),
      strip.text.y = element_text(angle = 0, color = "black", size = 13),
      strip.background = element_rect(fill = "gray90")
    )
  )
  
  # Manhattan Plot
  fused_fst <- ggplot() +
    # geom_point(data = fused_df, aes(x = midPos, y = Fst, color = centromere),
    #            alpha = 0.7, size = 0.9) +
    # scale_color_manual(values = c("maroon","gray20")) +
    geom_point(data = filter(fused_df,CategorizeSNP=='neither'), aes(x = midPos, y = Fst), 
               color = "maroon",alpha = 0.7, size = 0.9) +
    geom_point(data = filter(fused_df,CategorizeSNP=='fusion'), aes(x = midPos, y = Fst), 
               color = "gray30",alpha = 0.8, size = 1) +
    geom_point(data = filter(fused_df,CategorizeSNP=='centromere'), aes(x = midPos, y = Fst), 
               color = "black",alpha = 0.9, size = 1) +
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
    ylab(expression(italic(F[ST]))) +
    xlab("Position (Mb)") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(breaks = seq(10, 50, by = 10),
                       expand = expansion(mult = c(0.02,0.02)))
  #fused_fst
  
  # Save the plot to a pdf file
  jpeg(paste0("./figures/fst/",DATE,"/Pond-Chukchi_fused_chroms_fst_",format(Sys.Date(), "%Y%m%d"),".jpg"),
       width = 22, height = 10, res = 150, units = "in")
  print(fused_fst)
  dev.off()
  
  # Manhattan Plot
  fused_pi <- ggplot() +
    # geom_point(data = pi_df, aes(x = midPos, y = delta_pi, color = fusionPoint),
    #            alpha = 0.7, size = 1) +
    # scale_color_manual(values = c("mediumpurple3","gray20")) +
    geom_point(data = filter(pi_df,CategorizeSNP=='neither'), aes(x = midPos, y = delta_pi), color = "mediumpurple3",
               alpha = 0.5, size = 1.3,pch=16) +
    geom_point(data = filter(pi_df,CategorizeSNP=='fusion'), aes(x = midPos, y = delta_pi), color = "gray30",
               alpha = 0.6, size = 1.4,pch=19) +
    geom_point(data = filter(pi_df,CategorizeSNP=='centromere'), aes(x = midPos, y = delta_pi), color = "black",
               alpha = 0.7, size = 1.5,pch=19) +
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
    ylab(expression(Delta~pi)) +
    xlab("Position (Mb)") +
    scale_y_continuous(breaks = seq(-0.4, 0.4, by = 0.2),
                       limits = c(-0.48, 0.48),
                       expand = expansion(mult = c(0.01,0.01))) +
    scale_x_continuous(breaks = seq(10, 50, by = 10),
                       expand = expansion(mult = c(0.02,0.02)))
  fused_pi
  
  fused_fst <- fused_fst + theme(axis.title.x = element_blank(),
                                 axis.text.x = element_blank())
  
  # Save the plot to a pdf file
  jpeg(paste0("./figures/pi/boreogadus_Pond-Chukchi_fused_chroms_fst_pi_",format(Sys.Date(), "%Y%m%d"),".jpg"),
       width = 22, height = 10, res = 300, units = "in")
  print(fused_fst / fused_pi)
  dev.off()

  
### Only Chrs 1-3 #######################################
  
  fst3 <- ggplot() +
    geom_point(data = filter(fused_df, chr < 4), aes(x = midPos, y = Fst, color = centromere),
               alpha = 0.7, size = 0.9) +
    scale_color_manual(values = c("maroon","gray20")) +
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
    ylab(expression(italic(F[ST]))) +
    xlab("Position (Mb)") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_x_continuous(breaks = seq(10, 50, by = 10),
                       expand = expansion(mult = c(0.02,0.02)))+ 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  pi3 <- ggplot() +
    geom_point(data = filter(pi_df, chr < 4), aes(x = MB, y = delta_pi, color = centromere),
               alpha = 0.7, size = 1) +
    scale_color_manual(values = c("mediumpurple3","gray20")) +
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
    ylab(expression(Delta~pi)) +
    xlab("Position (Mb)") +
    scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.2),
                       limits = c(-0.5, 0.5)) +
    scale_x_continuous(breaks = seq(10, 50, by = 10),
                       expand = expansion(mult = c(0.02,0.02)))

  # Save the plot to a pdf file
  jpeg(paste0("./figures/pi/boreogadus_Pond-Chukchi_chroms1-3_fst_pi_",format(Sys.Date(), "%Y%m%d"),".jpg"),
       width = 22, height = 10, res = 150, units = "in")
  print(fst3 / pi3)
  dev.off()

  