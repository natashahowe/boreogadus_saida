# -- LD heatmap -------------------------------

# Code taken and edited from the [physalia tutorial](https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial3_ld_popstructure/scripts/LD_blocks.sh).

# --- Get LD HEATMAP from BioConductor --------------------
# Step 1, download BiocManager
# BiocManager::install(version = "3.18")
library(BiocManager)
# BiocManager::install(c("snpStats", "rtracklayer", "GenomicRanges", "GenomeInfoDb", "IRanges"))

# Install LDheatmap from github b/c it's no longer an r package
# install.packages("remotes")
# library(remotes)
# remotes::install_github("SFUStatgen/LDheatmap")

# Now can load LDheatmap
library(LDheatmap)

# --- After LDheatmap Install -------------------
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(scales)
library(here)
library(patchwork) 
library(LDheatmap)
library(reshape2)
library(gtools)

# --- Input Plink LD -------------------

CHROM="OZ177904.1"
PREFIX <- "boreogadus_OZ177904.1_plink"

# I had cut unnecessary columns from this ld output, so this won't be exactly the same code if it isn't edited
r <- read.table(paste0("./results/ld/", PREFIX,"_r2.ld"), header=F, stringsAsFactors=FALSE, 
                col.names = c('snp1', 'snp2', 'r2'))

id <- unique(mixedsort(c(r[,"snp1"],r[,"snp2"])))
posStart <- head(id,1)
posEnd <- tail(id,1)

# columns are snp1, snp2, and r2
# r <- rbind(r, c(posStart,posStart,0,NA,NA,NA,NA), c(posEnd,posEnd,0,NA,NA,NA,NA))
r <- rbind(r, c(posStart,posStart,NA), c(posEnd,posEnd,NA))

ld="r2"
m <- apply(acast(r, snp1 ~ snp2, value.var=ld, drop=FALSE),2,as.numeric)
rownames(m) <- colnames(m)
m <- m[mixedorder(rownames(m)),mixedorder(colnames(m))]
id <- rownames(m)
# I already removed the chromosome, so I just need to keep the dist value
# dist <- as.numeric(sub(".*:","",id))
dist <- as.numeric(id)
head(dist)
# Save plot
# Note: For dev.off() to work if running using ctrl+Enter, need to highlight pdf to dev.off() and run as one
pdf(paste0("../figures/melanostictis/blackspot-az_chr14_",posStart,"-",posEnd,"_LDheatmap_r2.pdf"), 
    width=10, height=10)
# png(paste0("../figures/melanostictis/blackspot-az_chr14_",posStart,"-",posEnd,"_LDheatmap_r2.jpeg"), 
#            width=900, height=900)
LDheatmap(m, genetic.distances=dist, 
          geneMapLabelX=0.5, geneMapLabelY=0.3, geneMapLocation = 0.15,
          title = "", color="blueToRed", LDmeasure="r",
          #name = "ld_plot"
)
dev.off()
