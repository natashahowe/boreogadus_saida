#######################################
### SCRIPT FOR CHECKING SYTENY OF REFERENCE GENOMES
### Sara Michele Schaal
### December 16, 2022
######################################


######################################
### INSTALL PACKAGES & LOAD FUNCTIONS
packages_needed <- c("ggplot2", "plotly", "ggpubr", "tidyverse","rsvg", "RIdeogram")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

g_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
######################################


######################################
## DIRECTORIES AND FILE NAMES

DATADIR <- "C:/Users/sara.schaal/Work/Pollock/data/"
SYNTENYDATA <- "genome/satsuma_summary.chained.out"
KARYOTYPEDATA <- "genome/karyotype_cod_pollock_synteny.csv"

######################################

df_synteny <- read.delim(paste0(DATADIR, SYNTENYDATA), sep = "\t", header = FALSE)
df_karyotype <- read.delim(paste0(DATADIR, KARYOTYPEDATA), sep = ",")
colnames(df_karyotype) <- c("Chr", "Start", "End", "fill", "species", "size", "color")

df_synteny_chroms_pcod <- df_synteny[grep("NC_", df_synteny$V1),]
df_synteny_chromsOnly <- df_synteny_chroms_pcod[complete.cases(df_synteny_chroms_pcod[grep("CM", df_synteny$V4),]),][,1:6]
colnames(df_synteny_chromsOnly) <- c("species1_chrom", "species1_start", "species1_end", "species2_chrom", 
                                     "species2_start", "species2_end")
df_synteny_chromsOnly$fill <- "cccccc"
df_synteny_chromsOnly$species1_chrom_num <- NULL
df_synteny_chromsOnly$species2_chrom_num <- NULL

for(i in 1:nrow(df_synteny_chromsOnly)){
  df_synteny_chromsOnly$species1_chrom_num[i] <-  str_extract(str_extract(df_synteny_chromsOnly$species1_chrom[i], "chromosome_\\d+"), "\\d+")
  df_synteny_chromsOnly$species2_chrom_num[i] <-  str_extract(str_extract(df_synteny_chromsOnly$species2_chrom[i], "chromosome_\\d+"), "\\d+")
}

#Species_1  Start_1    End_1 Species_2 Start_2   End_2   fill
df_synteny_order <- df_synteny_chromsOnly[c(8,2,3,9,5,6,7)]
colnames(df_synteny_order) <- c("Species_1", "Start_1", "End_1", "Species_2", "Start_2",   "End_2",   "fill")
df_synteny_order$Species_1 <- as.numeric(df_synteny_order$Species_1)
df_synteny_order$Species_2 <- as.numeric(df_synteny_order$Species_2)
df_synteny_order <- df_synteny_order[complete.cases(df_synteny_order),]
ideogram(karyotype = df_karyotype, synteny = df_synteny_order[seq(from = 1, to = nrow(df_synteny_order), by = 20),])
convertSVG("chromosome.svg", device = "png")

