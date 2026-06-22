# Pi for Each Pop


### INSTALL PACKAGES & LOAD FUNCTIONS ######################################

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "tools", 
                     "here", "patchwork", "readxl","KernSmooth")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

DATE="20241222" 

### Read in tab-delimited table that has two columns chrom name from angsd and a simplified name (e.g., chr_1)
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

##### Pi #####################################################

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
  left_join(chrom_df)

# add fusionPoint locations
pi_df <- pi_df %>%
  rowwise() %>%
  mutate(fusionPoint = any(chr == fusionPoint$chr &
                            midPos >= fusionPoint$minPos &
                            midPos <= fusionPoint$maxPos),
         centromere = any(chr == centromere$chr &
                            midPos >= centromere$minPos &
                            midPos <= centromere$maxPos))

head(pi_df)

# priority goes to centromere
pi_df <- pi_df %>%
  mutate(CategorizeSNP = case_when(centromere == T ~ "centromere",
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
    axis.title.x = element_text(size = 20, color = "black"),
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


# Pi Pop1
fused_pi1 <- ggplot() +
  # geom_point(data = pi_df, aes(x = midPos, y = Pond_pi, color = CategorizeSNP),
  #            alpha = 0.7, size = 1.3) +
  #geom_rect(data = fusionPoint, aes(xmin=minPos,xmax=maxPos,ymin=-Inf,ymax=Inf),alpha=0.5,fill="gray")+
  #geom_rect(data = centromere, aes(xmin=minPos,xmax=maxPos,ymin=-Inf,ymax=Inf),alpha=0.5,fill="black")+
  geom_point(data = filter(pi_df,CategorizeSNP=='neither'), aes(x = midPos, y = Pond_pi), color = "darkorange3",
             alpha = 0.7, size = 1.5,pch=16) +
  geom_point(data = filter(pi_df,CategorizeSNP=='fusion'), aes(x = midPos, y = Pond_pi), color = "gray25",
             alpha = 0.7, size = 1.6,pch=19) +
  geom_point(data = filter(pi_df,CategorizeSNP=='centromere'), aes(x = midPos, y = Pond_pi), color = "black",
             alpha = 0.9, size = 1.8,pch=19) +
  #scale_color_manual(values = c("gray40","black","darkorange3")) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  ylab(expression(Pond~(pi))) +
  xlab("Position (Mb)") +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1),
                     limits = c(0, 0.52)) +
  scale_x_continuous(breaks = seq(10, 50, by = 10),
                     expand = expansion(mult = c(0.02,0.02)))
# fused_pi1

# Pi Pop2
fused_pi2 <- ggplot() +
  geom_point(data = filter(pi_df,CategorizeSNP=='neither'), aes(x = midPos, y = Chukchi_pi), color = "maroon4",
             alpha = 0.7, size = 1.5,pch=16) +
  geom_point(data = filter(pi_df,CategorizeSNP=='fusion'), aes(x = midPos, y = Chukchi_pi), color = "gray25",
             alpha = 0.7, size = 1.6,pch=19) +
  geom_point(data = filter(pi_df,CategorizeSNP=='centromere'), aes(x = midPos, y = Chukchi_pi), color = "black",
             alpha = 0.9, size = 1.8,pch=19) +
  #scale_color_manual(values = c("maroon4","black")) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  ylab(expression(Chukchi~(pi))) +
  xlab("Position (Mb)") +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1),
                     limits = c(0, 0.52)) +
  scale_x_continuous(breaks = seq(10, 50, by = 10),
                     expand = expansion(mult = c(0.02,0.02)))
# fused_pi2


fused_pi <- (fused_pi1 + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank())) / 
  (fused_pi2)
# fused_pi

# Save the plot to a pdf file
jpeg(paste0("./figures/pi/boreogadus_Pond-Chukchi_fused_chroms_pi_",format(Sys.Date(), "%Y%m%d"),".jpg"),
     width = 18, height = 9, res = 300, units = "in")
print(fused_pi)
dev.off()

### EXTRA DELTA PI ADD ON #################################################

# Delta Pi
delta_pi <- ggplot() +
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
# delta_pi


fused_delta_pi <- (fused_pi1 + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank())) / 
  (fused_pi2 + theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank())) / 
  delta_pi

# Save the plot to a pdf file
jpeg(paste0("./figures/pi/boreogadus_Pond-Chukchi_fused_chroms_pi_",format(Sys.Date(), "%Y%m%d"),".jpg"),
     width = 18, height = 12, res = 300, units = "in")
print(fused_delta_pi)
dev.off()

#### TEST OUT KERNEL SMOOTHING #################################################
#install.packages('smoothr')
library(smoothr)

Chukchi_pi <- pi_df %>%
  filter(!is.na(Chukchi_pi)) %>%  # have to filter out NA for smooth_ksmooth() fxn to work
  select(chr, midPos, Chukchi_pi) %>%
  rename(pi = Chukchi_pi)

# apply the kernel smoothing
Chukchi_smooth <- Chukchi_pi %>%
  as_tibble() %>%
  group_split(chr) %>%
  map_dfr(~ {
    mat <- as.matrix(select(.x, midPos, pi))
    smoothed <- smooth_ksmooth(mat, wrap = TRUE, bandwidth = 2)
    as.data.frame(smoothed) %>%
      rename(midPos = V1, pi = V2) %>%
      mutate(chr = unique(.x$chr))
  })


class(Chukchi_pi)
class(Chukchi_smooth)
#plot(Chukchi_smooth, type = "l", col = "black", lwd = 1)
#lines(m_smooth, lwd = 3, col = "red")

ggplot() +
  geom_line(data = as.data.frame(Chukchi_smooth), aes(x = midPos, y = pi), color='black') +
  facet_wrap(. ~chr, scales = "free") +
  labs(y = 'Pi',
       x = "Chromosome Position (Mb)") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     #limits = c(xstart, xend),
                     breaks = seq(0, 100, by = 10))


Pond_pi <- pi_df %>%
  filter(!is.na(Pond_pi)) %>%  # have to filter out NA for smooth_ksmooth() fxn to work
  select(chr, midPos, Pond_pi) %>%
  rename(pi = Pond_pi)

# apply the kernel smoothing
Pond_smooth <- Pond_pi %>%
  as_tibble() %>%
  group_split(chr) %>%
  map_dfr(~ {
    mat <- as.matrix(select(.x, midPos, pi))
    smoothed <- smooth_ksmooth(mat, wrap = TRUE, bandwidth = 2)
    as.data.frame(smoothed) %>%
      rename(midPos = V1, pi = V2) %>%
      mutate(chr = unique(.x$chr))
  })

# Pi Pop1
fused_pi1 <- ggplot() +
  geom_point(data = pi_df, aes(x = midPos, y = Pond_pi, color = fusionPoint),
             alpha = 0.7, size = 1.3) +
  scale_color_manual(values = c("darkorange3","gray40")) +
  geom_line(data = Pond_smooth, aes(x = midPos, y = pi), 
            linewidth = 0.90, color='black') +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  ylab(expression(Pond~(pi))) +
  xlab("Position (Mb)") +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1),
                     limits = c(0, 0.52)) +
  scale_x_continuous(breaks = seq(10, 50, by = 10),
                     expand = expansion(mult = c(0.02,0.02)))
fused_pi1

# Pi Pop2
fused_pi2 <- ggplot() +
  geom_point(data = pi_df, aes(x = midPos, y = Chukchi_pi, color = fusionPoint),
             alpha = 0.7, size = 1.3) +
  scale_color_manual(values = c("maroon","gray40")) +
  geom_line(data = Chukchi_smooth, aes(x = midPos, y = pi), 
            linewidth = 0.90, color='black') +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  ylab(expression(Chukchi~(pi))) +
  xlab("Position (Mb)") +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1),
                     limits = c(0, 0.52)) +
  scale_x_continuous(breaks = seq(10, 50, by = 10),
                     expand = expansion(mult = c(0.02,0.02)))
#fused_pi2


fused_pi_smooth <- (fused_pi1 + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank())) / 
  (fused_pi2)
fused_pi_smooth

# Save the plot to a pdf file
jpeg(paste0("./figures/pi/boreogadus_Pond-Chukchi_fused_chroms_pi_smoothed_line_",format(Sys.Date(), "%Y%m%d"),".jpg"),
     width = 18, height = 9, res = 150, units = "in")
print(fused_pi_smooth)
dev.off()

### Tajima D ###################################################


tajima <- read_xlsx("./results/pi/Tajimas_D_Pond_Chukchi.xlsx") %>%
  rename(chrName = Chr) %>%
  left_join(chrom_df, by = "chrName") %>%
  relocate(chr, .before=chrName)

# Plot points, Tajima
tajima %>%
  pivot_longer(cols = c(Tajima_Chukchi,Tajima_Pond), 
               names_to = "Pop", values_to = "tajD") %>%
  mutate(Pop=sub("Tajima_","",Pop),
         midPos = WinCenter/1e6) %>%
  filter(chr < 6) %>%
  ggplot() +
  geom_point(aes(x = midPos, y = tajD), color = "gray30",
             alpha = 0.7, size = 1.3) +
  facet_grid(Pop ~ chr, scales = "free_x", space = "free_x") +
  ylab("Tajima's D") + xlab("Position (Mb)")


######  Kernel - BW5   ##################
library(KernSmooth)

# apply the kernel smoothing
Taj_kernel <- tajima %>%
  pivot_longer(cols = c(Tajima_Chukchi,Tajima_Pond), 
               names_to = "Pop", values_to = "tajD") %>%
  mutate(Pop=sub("Tajima_","",Pop),
         midPos = WinCenter/1e6) %>%
  filter(chr < 6)%>%
  drop_na(chr,midPos, tajD) %>%
  dplyr::select(chr, midPos, tajD,Pop) %>%
  as_tibble() %>%
  group_by(chr,Pop) %>%
  group_modify(~ {
    if (nrow(.x) < 2) return(tibble(midPos = numeric(), Taj_Kdensity = numeric()))
    
    x_data <- .x$midPos
    y_data <- .x$tajD
    
    smoothed <- locpoly(
      x = x_data, y = y_data,
      bandwidth = 5,
      degree = 0,
      range.x = range(x_data),
      gridsize = 4001
    )
    
    tibble(
      midPos = smoothed$x,
      Taj_Kdensity = smoothed$y
    )
  }) %>%
  ungroup()

# Plot kernel lines
tajd <- ggplot() +
  geom_line(data = as.data.frame(Taj_kernel), aes(x = midPos, y = Taj_Kdensity, color=Pop),
            size=1.5)+
  facet_grid(. ~chr, scales = "free_x", space = "free_x") +
  scale_color_manual(name="Population",values=c("maroon4","darkorange3"))+
  labs(y = "Tajima's D",
       x = "Chromosome Position (Mb)") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(breaks = c(-0.4, -0.2,0, 0.2,0.4,0.6,0.8))+
  theme(legend.position = c(0.06,0.9),
        legend.text = element_text(size=16))


# Save the plot to a pdf file
jpeg(paste0("./figures/pi/boreogadus_Pond-Chukchi_fused_chroms_tajD_kernel_bw5_",format(Sys.Date(), "%Y%m%d"),".jpg"),
     width = 18, height = 5, res = 150, units = "in")
print(tajd)
dev.off()

### Combine Pi and Taj D Plots ####################

# Add Tajima's D to the end of the pop=specific pi
pi_plus_taj <- (fused_pi1 + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank())) / 
  (fused_pi2 + theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank())) / 
  (tajd)
pi_plus_taj

# Save the plot to a pdf file
jpeg(paste0("./figures/pi/boreogadus_Pond-Chukchi_fused_chroms_pi_plus_tajD_kernel_bw5_",format(Sys.Date(), "%Y%m%d"),".jpg"),
     width = 18, height = 16, res = 150, units = "in")
print(pi_plus_taj)
dev.off()
