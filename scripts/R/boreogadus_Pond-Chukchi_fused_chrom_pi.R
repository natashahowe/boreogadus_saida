# Pi for Each Pop

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
# Read in tab-delimited table that has two columns chrom name from angsd and a simplified name (e.g., chr_1)
chrom_df <- read.table("./data/R/Arctic_cod_genome_chromosomes.txt", header = T)

# Read in metadata file with population, sampling location, and region
meta_df <- read.delim("./data/R/region_metadata.txt", header = T)

# read in color file
color_df <- read.delim("./data/R/color_metadata_allpops_flip.txt", header = T)

centromere <- data.frame(chr = 1:5,
                         minPos = c(24,28,36,22,18.5),
                         maxPos = c(30,37,40,24,20))

##### Pi #####################################################

PIFILES <- Sys.glob(file.path(here(),"results","pi","pop_pi","b*")) 
pi_dfs <- NULL
for(i in 1:length(PIFILES)){
  pi_dfs[[i]] <- read_csv(PIFILES[i], col_types = cols())
}

head(pi_dfs[[1]])

do.call(rbind,pi_dfs) %>%
  group_by(chrName) %>%
  summarize(n())

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

# add centromere locations
pi_df <- pi_df %>%
  rowwise() %>%
  mutate(centromere = any(chr == centromere$chr &
                            midPos >= centromere$minPos &
                            midPos <= centromere$maxPos))

head(pi_df)

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


# Pi Pop1
fused_pi1 <- ggplot() +
  geom_point(data = pi_df, aes(x = midPos, y = Pond_pi, color = centromere),
             alpha = 0.7, size = 1.3) +
  scale_color_manual(values = c("darkorange3","black")) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  ylab(expression(Pond~(pi))) +
  xlab("Position (Mb)") +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1),
                     limits = c(0, 0.52)) +
  scale_x_continuous(breaks = seq(10, 50, by = 10),
                     expand = expansion(mult = c(0.02,0.02)))
#fused_pi1

# Pi Pop2
fused_pi2 <- ggplot() +
  geom_point(data = pi_df, aes(x = midPos, y = Chukchi_pi, color = centromere),
             alpha = 0.7, size = 1.3) +
  scale_color_manual(values = c("maroon4","black")) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  ylab(expression(Chukchi~(pi))) +
  xlab("Position (Mb)") +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1),
                     limits = c(0, 0.52)) +
  scale_x_continuous(breaks = seq(10, 50, by = 10),
                     expand = expansion(mult = c(0.02,0.02)))
#fused_pi2


fused_pi <- (fused_pi1 + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank())) / 
  (fused_pi2)
fused_pi

# Save the plot to a pdf file
jpeg(paste0("./figures/pi/boreogadus_Pond-Chukchi_fused_chroms_pi_",format(Sys.Date(), "%Y%m%d"),".jpg"),
     width = 18, height = 9, res = 150, units = "in")
print(fused_pi)
dev.off()

### EXTRA DELTA PI ADD ON #################################################

# Delta Pi
delta_pi <- ggplot() +
  geom_point(data = pi_df, aes(x = midPos, y = delta_pi, color = centromere),
             alpha = 0.7, size = 1) +
  scale_color_manual(values = c("mediumpurple3","gray20")) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  ylab(expression(Delta~pi)) +
  xlab("Position (Mb)") +
  scale_y_continuous(breaks = seq(-0.4, 0.4, by = 0.2),
                     limits = c(-0.48, 0.48),
                     expand = expansion(mult = c(0.01,0.01))) +
  scale_x_continuous(breaks = seq(10, 50, by = 10),
                     expand = expansion(mult = c(0.02,0.02)))
delta_pi


fused_delta_pi <- (fused_pi1 + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank())) / 
  (fused_pi2 + theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank())) / 
  delta_pi

# Save the plot to a pdf file
jpeg(paste0("./figures/pi/boreogadus_Pond-Chukchi_fused_chroms_pi_",format(Sys.Date(), "%Y%m%d"),".jpg"),
     width = 18, height = 12, res = 150, units = "in")
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
  #geom_line(data = chum_smooth, aes(x = midPos, y = Fst), color = chumcol) +
  #geom_line(data = sock_smooth, aes(x = midPos, y = Fst), color = sockcol) +
  facet_wrap(. ~chr, scales = "free") +
  labs(y = 'Pi',
       x = "Chromosome Position (Mb)") +
  #ggtitle("Lrrc9 Gene Region") +
  # scale_y_continuous(#limits = c(-0.02, 1.02),
  #                    breaks = seq(-10, 10, by = 1),
  #                    expand = expansion(mult = c(0.001, 0.01))) +
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
  geom_point(data = pi_df, aes(x = midPos, y = Pond_pi, color = centromere),
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
  geom_point(data = pi_df, aes(x = midPos, y = Chukchi_pi, color = centromere),
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

