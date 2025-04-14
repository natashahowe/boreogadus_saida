######################################
### INSTALL PACKAGES & LOAD FUNCTIONS

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "tools", "here")

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
  # Read in tab-delimited table that has two columns chrom name from angsd and a simplified name (e.g., chr_1)
  chrom_df <- read.table("./data/R/Arctic_cod_genome_chromosomes.txt", header = T)

  # Read in metadata file with population, sampling location, and region
  meta_df <- read.delim("./data/R/region_metadata.txt", header = T) #%>%
    #add_row(Name = "Chukchi Subset", Pop = "Chukchi_subset", Region = "West Arctic") %>%
    #add_row(Name = "Pond Subset", Pop = "Pond_subset", Region = "East Arctic")
  
  # read in color file
  color_df <- read.delim("./data/R/color_metadata_allpops_flip.txt", header = T) #%>%
    #add_row(Population = "Chukchi Subset", Pop = "Chukchi_subset", Region = "West Arctic", Color = "maroon", Color2 = "maroon") %>%
    #add_row(Population = "Pond Subset", Pop = "Pond_subset", Region = "East Arctic", Color = "darkorange3", Color2 = "darkorange3")
  
  # Specify the order of some factors for plotting later
  meta_df$Pop <- factor(meta_df$Pop, levels = meta_df$Pop)
  meta_df$Region <- factor(meta_df$Region, levels = unique(meta_df$Region))
  
  POPLIST <- unique(meta_df$Pop)
  POPLIST
  
  # create mypalette
  mypalette <- as.vector(color_df$Color) # turn colors into vector
    names(mypalette) <- color_df$Pop # attach pop name to palette color
    mypalette
  
  
  POPLIST <- POPLIST[-grep('Labrador')]  
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
    select(chr, midPos, Fst, comparison) %>%
    mutate(midPos = as.numeric(midPos)/1e6,
           Fst = as.numeric(Fst))
  head(fst_df)
  
  # Strip the original file name such that only the population that FOCALPOP is 
    # being compared to is retained (for plot labeling)
  final_df <- fst_df %>%
    mutate(temp = basename(tools::file_path_sans_ext(comparison)),
           temp = sub("boreogadus_","",temp),
           temp = sub("_subset","",temp),
           temp = sub("_downsampled","",temp),
           temp = gsub(paste0("_",DATE,".fst.SNP"), "", temp),
           temp = gsub(FOCALPOP, "", temp),
           POP2 = gsub("-", "", temp)) %>%
    select(-c(comparison, temp))
  
  head(final_df)
  unique(final_df$POP2)
  
    rm(fst_df)
  
  # Make the POP2 and chr column in the data frame into a factor with a specific order
  final_df$POP2 <- factor(final_df$POP2, levels = unique(meta_df$Pop))
  final_df$chr <- factor(final_df$chr, levels = chrom_df$chr)
  
    # get rid of negative fst values
    final_df$Fst[final_df$Fst < 0] <- 0
  
  #################################################################
  # Part 2: Plot the data
  
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
  
  # Manhattan Plot
  manhplot <- ggplot() +
    geom_point(data = final_df, aes(x = midPos, y = Fst, color = POP2), 
               alpha = 0.9, size = 0.9) +
    scale_color_manual(values = mypalette) +
    facet_grid(POP2 ~ chr, scales = "free_x", space = "free_x",
               labeller = labeller(facet_category = label_wrap_gen(width = 13))) +
    ylab(expression(italic(F[ST]))) +
    xlab("Position (Mb)") +
    ggtitle(TITLE) +
    scale_y_continuous(breaks = seq(0.25, 1, by = 0.25)) +
    scale_x_continuous(breaks = seq(10,40, by = 10))
  
  # Save the plot to a pdf file
  jpeg(paste0("./figures/fst/",DATE,"/",FOCALPOP, "_fst_wholegenome_noLabrador.jpg"),
       width = 22, height = 10, res = 150, units = "in")
  print(manhplot)
  dev.off()
  
  rm(manhplot)

  for(j in 1:nrow(chrom_df)){
    chrom = chrom_df$chr[j]

    # removed chr level plots for now
    chr_df <- final_df %>% 
      filter(chr == chrom)
    
    # Manhattan Plot
    chrplot <- ggplot() +
      geom_point(data = chr_df, aes(x = midPos, y = Fst, color = POP2), 
                 alpha = 0.9, size = 1.5) +
      scale_color_manual(values = mypalette) +
      facet_grid(POP2 ~ ., scales = "free_x", space = "free_x",
                 labeller = labeller(facet_category = label_wrap_gen(width = 13))) +
      ylab(expression(italic(F[ST]))) +
      xlab("Position (Mb)") +
      ggtitle(paste0(TITLE,": Chromosome ",chrom)) +
      scale_y_continuous(breaks = seq(0.25, 1, by = 0.25)) +
      scale_x_continuous(breaks = seq(0, 60, by = 5),
                         expand = c(0.01,0.001,0.001,0.01))
    
    # Save the plot to a pdf file
    jpeg(paste0("./figures/fst/",DATE,"/",FOCALPOP,"/",FOCALPOP,"_noLabrador_fst_chr",chrom,".jpg"),
         width = 18, height = 12, res = 150, units = "in")
    print(chrplot)
    dev.off()
  } 
  rm(final_df, chrplot)
}  
  # END OF GENOME PLOT LOOP
 
#####################################################################
# chr BY chr
#####################################################################


  # Set the chr theme
  theme_set(
    theme( 
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(angle = 0, size = 12, vjust = 0.5, color = "black"),
      axis.title = element_text(size = 20, color = "black"),
      title = element_text(size = 16, color = "black"),
      panel.background = element_rect(fill = "white"), 
      panel.spacing = unit(0,"lines"),
      strip.text.x = element_text(angle = 0, color = "black", size = 14),
      strip.text.y = element_text(angle = 0, color = "black", size = 13)
    )
  )
  
  for(j in 1:length(unique(final_df$chr))){
    plotchr <- levels(final_df$chr)[j]
    TITLE <- meta_df$Name[which(meta_df$Pop == FOCALPOP)]
    chrplot <- ggplot() +
      geom_point(data = filter(final_df, chr == plotchr), aes(x = midPos, y = Fst, color = POP2), 
                 alpha = 1, size = 1.2) +
      scale_color_manual(values = mypalette) +
      facet_grid(POP2 ~ ., scales = "free_x") +
      ylab(expression(italic(F[ST]))) +
      xlab("Position (Mb)") +
      ggtitle(paste0(TITLE,": chr",plotchr)) +
      # set tick mark spacing
      scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                         limits = c(0, 0.82)) +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                         breaks = seq(0, 50, by = 5))
    ggsave(paste0("./figures/fst/chr/",FOCALPOP,"/",FOCALPOP,"_chr",plotchr, "_fst.jpg"), 
           plot = chrplot, width = 15, height = 15, units = "in")
  }
  