# plot pink, sockeye, chum aligned to chum along with genes in the region
# kernel smoothing
# 06/05/2024
# chr 35 at lrrc9

packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "tidyverse",
                     "stringr", "data.table", "plyr","gtools","reshape2", 
                     "patchwork", "RColorBrewer", "smoothr")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}
########################################################

# which region of chr35 to plot (first and last position)
# panel spanning larger region (Fig 5)
xstart = 27.86
xend = 28.3

# lrrc9 only - Zoomed in panel (Suppl. Fig)
#xstart = 28.127
#xend = 28.17

############ IMPORT GENE DATA FROM NCBI ###################################

# find exons from gff file for genes of interest (from NCBI chum reference genome)
gff_df <- read.delim("./data/R/genomic.gff", header = F, comment.char = "#")

# remove excess columns
gff_df <- gff_df[,c(1:5,9)]

# rename remaining columns
colnames(gff_df) <- c("chrName", "RefSeq","exon","start.pos","fin.pos", "ID")

# only keep chr35
gff_chr35 <- gff_df %>%
  filter(chrName == "NC_068455.1")
rm(gff_df)

############ IMPORT AND EDIT SPECIES ALLELE FST FILES #####################  

########## PINK
pink_Fst <- read.delim2("./results/fst/allele/pink-chum_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
colnames(pink_Fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

# create chr column
pink_Fst$chr <- 35

pink_Fst <- pink_Fst %>%
  mutate(midPos = midPos/1e6) %>%
  select(chr, midPos, Fst)

# remove negative Fst values
pink_Fst$Fst[pink_Fst$Fst < 0] <- 0
head(pink_Fst)

# Cut Fst to start and end
pink_df <- pink_Fst %>%
  filter(midPos > xstart - 0.5,
         midPos < xend + 0.5)

pink_df$Fst <- as.numeric(pink_df$Fst)
pink_df$midPos <- as.numeric(pink_df$midPos)

############### SOCKEYE
sock_fst <- read.delim("./results/fst/allele/sock-chum_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                       row.names = NULL,sep = "\t")
colnames(sock_fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

# create chr column
sock_fst$chr <- 35

# change midPos to Mb, and remove unnecessary columns
sock_fst <- sock_fst %>%
  mutate(midPos = midPos/1e6) %>%
  select(chr, midPos, Fst)

# remove negative Fst values
sock_fst$Fst[sock_fst$Fst < 0] <- 0
head(sock_fst)

sock_df <- sock_fst %>%
  filter(midPos > xstart - 0.5,
         midPos < xend + 0.5)

sock_df$Fst <- as.numeric(sock_df$Fst)
sock_df$midPos <- as.numeric(sock_df$midPos)

########## CHUM
chum_Fst <- read.delim2("./results/fst/allele/chumrun_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
colnames(chum_Fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

# create chr column
chum_Fst$chr <- 35

chum_Fst <- chum_Fst %>%
  mutate(midPos = midPos/1e6) %>%
  select(chr, midPos, Fst)

# remove negative Fst values
chum_Fst$Fst[chum_Fst$Fst < 0] <- 0
head(chum_Fst)

# filter both pops to desired start and end point
chum_df <- chum_Fst %>%
  filter(midPos > xstart - 0.5,
         midPos < xend + 0.5) 

chum_df$Fst <- as.numeric(chum_df$Fst)
chum_df$midPos <- as.numeric(chum_df$midPos)

############ ADD GENES #################################

# only the region of interest
gff_region <- gff_chr35 %>%
  mutate(start.pos = start.pos/1e6,
         fin.pos = fin.pos/1e6) %>% 
  filter(fin.pos > xstart,  
         start.pos < xend)

# prep pattern for str_match below
gene_pattern <- "gene=\\s*(.*?)\\s*;"    # keep string btwn "gene=" & ":product" 
exon_pattern <- "ID=exon-\\s*(.*?)\\s*;Parent"  # keep string btwn "exon=" & ";Parent" 
descr_pattern <- ";description=\\s*(.*?)\\s*;"

# for first creating gene file
gene_write35 <- gff_region %>% 
  filter(exon == "gene") %>%
  mutate(gene = str_match(ID, gene_pattern)[,2],
         geneName = str_match(ID, descr_pattern)[,2])

#write.csv(gene_write35, "./data/R/genomic_chr35_lrrc9_region.csv", row.names = F)

# create new columns for genes and exons from ID
gff_region_exon <- gff_region %>%
  filter(exon == "exon") %>%
  mutate(gene = str_match(ID, gene_pattern)[,2],  # gene abbr.      
         exonID = str_match(ID,exon_pattern)[,2]) # mRNA name and exon number

# this exon in six6a is too small that it doesn't even plot
# make slightly larger so it is visible in plot
gff_region_exon$fin.pos[which(gff_region_exon$gene == "six6a")[1]] <- 27.994500 # changed from 27.994277

# only retain columns of interest
exons_to_plot <- gff_region_exon[,c(4,5,7:8)]
unique(exons_to_plot$gene)

# factor based on gene name
exons_to_plot$gene <- factor(exons_to_plot$gene, levels = unique(exons_to_plot$gene))
levels(exons_to_plot$gene)

########### GENE & COLORS - ACTS AS INTRONS ############
# this file is manually edited for plotting purposes
# colors are assigned to each gene, two versions for different panel sizes
if(xstart == 28.127){ # zoomed in lrrc9 plot
  genes_df <- read.delim2("./data/R/chum_genes_exons2.txt", header = T, 
                          sep = "\t", row.names = NULL)
}else{ # larger plot
  genes_df <- read.delim2("./data/R/chum_genes_exons.txt", header = T, 
                          sep = "\t", row.names = NULL) 
}

# set factors for plotting columns
genes_df$gene <- factor(genes_df$gene, levels = genes_df$gene)

# assign name of gene to the color so it can be plotted properly
mypalette <- genes_df$color
names(mypalette) <- levels(genes_df$gene)
mypalette

# having an issue with these not being numeric, so assign them all as such here
genes_df$beg.pos = as.numeric(genes_df$beg.pos)
genes_df$end.pos = as.numeric(genes_df$end.pos)
genes_df$y.min = as.numeric(genes_df$y.min)
genes_df$y.max = as.numeric(genes_df$y.max)

# only keep exons from genes that have color codes
# this removes all that start with LOC####
exons_to_plot <- filter(exons_to_plot, gene %in% genes_df$gene)
head(exons_to_plot)

################# PLOTTING ###########################

# Species Specific Colors
spp_palette <- palette.colors(palette = "R4")[1:4]
spp_palette <- c("#000000", "#DF536B", "#61D04F", "#2297E6")
names(spp_palette) <- c("Coho", "Pink", "Chum", "Sockeye")
pinkcol <- palette.colors(palette = "R4")[2]
chumcol <- palette.colors(palette = "R4")[3]
sockcol <- palette.colors(palette = "R4")[4]
  
# Set the general themes
theme_set(
  theme( 
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor.x = element_line(color = "gray90"),
    axis.text.y = element_text(angle = 0, size = 11, color = "black", vjust = 0.5),
    axis.title.y = element_text(size = 13, angle = 90),
    strip.text.y = element_text(angle = 0),
    panel.background = element_rect(fill = "white"), 
    axis.line = element_line(),
    panel.border = element_rect(color = "black", fill = "NA")
  )
)


pink_lrrc9 <- as.matrix(pink_df[,c("midPos", "Fst")])
pink_smooth <- smooth_ksmooth(pink_lrrc9, wrap = TRUE, bandwidth = 3)
colnames(pink_smooth) <- c("midPos", "Fst")
class(pink_lrrc9)
class(pink_smooth)
#plot(pink_smooth, type = "l", col = "black", lwd = 1)
#lines(m_smooth, lwd = 3, col = "red")

chum_lrrc9 <- as.matrix(chum_df[,c("midPos", "Fst")])
chum_smooth <- smooth_ksmooth(chum_lrrc9, wrap = TRUE, bandwidth = 3)
colnames(chum_smooth) <- c("midPos", "Fst")

sock_lrrc9 <- as.matrix(sock_df[,c("midPos", "Fst")])
sock_smooth <- smooth_ksmooth(sock_lrrc9, wrap = TRUE, bandwidth = 3)
colnames(sock_smooth) <- c("midPos", "Fst")

pink_plot <- ggplot() +
  geom_line(data = pink_smooth, aes(x = midPos, y = Fst), color = pinkcol) +
  geom_line(data = chum_smooth, aes(x = midPos, y = Fst), color = chumcol) +
  geom_line(data = sock_smooth, aes(x = midPos, y = Fst), color = sockcol) +
  labs(y = expression(F[ST]),
       x = "Chromosome Position (Mb)") +
  ggtitle("Lrrc9 Gene Region") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.2),
                     expand = expansion(mult = c(0.001, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart, xend),
                     breaks = seq(0, 100, by = 0.1))
pink_plot

jpeg("./figures/fst/testkernel1.jpg", quality=100, 
     height=1000, width=3000, pointsize=14, res=600)
plot(chum_smooth, type = "l", col = chumcol, lwd = 2)
lines(pink_smooth, lwd = 2, col = pinkcol)
lines(sock_smooth, lwd = 2, col = sockcol)
dev.off()

################################################################
################################################################
sock_plot <- ggplot() +
  geom_point(data = sock_df, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.5, color = "gray10") + #01665E
  labs(y = expression(F[ST]~-~Sockeye)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(-0.02, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart, xend),
                     breaks = seq(0, 100, by = 0.1)) +           # change from 0.1 to 0.01 for zoom !!!
  theme(
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor.x = element_line(color = "gray90"),
    axis.text.y = element_text(angle = 0, size = 12, color = "black", vjust = 0.5),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 13, angle = 90),
    axis.title.x = element_blank(),
    strip.text.y = element_text(angle = 0),
    panel.background = element_rect(fill = "white"), 
    axis.line = element_line(),
    plot.margin = unit(c(0,0.15,0.02,0.05), "cm"),
    panel.border = element_rect(color = "black", fill = "NA"))

chum_plot <- ggplot() +
  geom_point(data = chum_df, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.5, color = "gray10") +  #5D478B
  labs(x="Chromosome Position (Mb)", 
       y = expression(F[ST]~-~Chum)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart, xend),
                     breaks = seq(0, 100, by = 0.1)) +
  theme(
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor.x = element_line(color = "gray90"),
    axis.text = element_text(angle = 0, size = 12, color = "black", vjust = 0.5),
    axis.title.y = element_text(size = 13, angle = 90),
    axis.title.x = element_text(size = 12),
    strip.text.y = element_text(angle = 0),
    panel.background = element_rect(fill = "white"), 
    axis.line = element_line(),
    plot.margin = unit(c(0,0.15,0.05,0.05), "cm"),
    panel.border = element_rect(color = "black", fill = "NA"))
#chum_plot

####### plot genes/exons
# optional add on arrows specifically for lrrc9 region (for zoomed plot)
#   uncomment geom_segment and scale_color_manual
gene_plot <- ggplot() +
  geom_rect(data = genes_df, aes(xmin = beg.pos, xmax = end.pos, 
                                 ymin = y.min, ymax = y.max,
                                 fill = gene)) +
  # uncomment next two lines to add arrows
  #  geom_segment(data = genes_df, aes(x = end.pos, xend = ((end.pos - beg.pos)/2 + beg.pos),y = (y.max - y.min)/2 + y.min, yend = (y.max - y.min)/2 + y.min,
  #                                  color = gene), arrow = arrow(type = "closed", length = unit(0.15, "inches")), show.legend = F) +
  geom_rect(data = exons_to_plot, aes(xmin = start.pos, xmax = fin.pos, 
                                      ymin = 0, ymax = 0.1,
                                      fill = gene)) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart, xend)) +
  scale_fill_manual(values = mypalette) +
  # uncomment for arrows
  # scale_color_manual(values = mypalette) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1,"lines"),
    axis.line = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.margin = unit(c(0.05,0,0,0), "cm")) +
  guides(fill = guide_legend(title = "Genes"))

# plot three figures on top of one another
multiplot <- pink_plot / gene_plot / sock_plot / gene_plot / chum_plot + 
  plot_layout(heights = c(1, 0.13, 1, 0.13, 1),
              guides = "collect")
multiplot

# output figure
# change output name if switching between zoomed in and non-zoomed figures !!!
jpeg(paste0("./figures/fst/threespp_lrrc9_",xstart,"-",xend,"Mb_fst_genes.jpg"), 
     width = 12, height = 8, res = 150, units = "in")
print(multiplot)
dev.off()

# END OF PLOT SCRIPT
#############################################################################