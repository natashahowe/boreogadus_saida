# PCA for Sex-Associated region delineated in Hoff

library(here)
library(tools)
library(tidyverse)
library(patchwork)
library(ggnewscale)

#### OPTIONS #########################################################
suffix = '_20241222'
WILDCARD = '*908.1_s*e*.cov' # sex PCA
COVFILES <- Sys.glob(file.path(here(),"results","pca",WILDCARD))
basename(COVFILES)

## Filepaths ################################################################
BAMFILE <- "./data/R/boreogadus_filtered_downsampled_bams.txt"
METADATAFILE <- "./data/R/boreogadus_metadata.csv"
COLORFILE <- "./data/R/color_metadata_allpops_flip.txt"
chrom_df <- read.table('./data/R/boreogadus_chromosomes.txt', header = T)

### Metadata ###############################################################
meta_df <- read.csv(METADATAFILE, header = T) %>%
  mutate(Population = ifelse(Population=='Southhampton Island', 'Southampton Island', Population),
         Pop = case_when(Population == 'SE Baffin' ~ 'SEBaffin',
                         TRUE ~ gsub(' .*$','',Population))) %>%
  filter(Pop != 'Unknown')

color_df <- read.delim2(COLORFILE, header = T, row.names = NULL, sep = "\t") %>%
  mutate(Pop = gsub('hh','h',Pop)) %>%
  left_join(meta_df %>% distinct(Population, Pop))

# call in bams
bam_df <- read.table(BAMFILE, header = F) %>%
  mutate(sampleID = basename(V1) %>% tools::file_path_sans_ext(.),
         sampleID = gsub("_sorted",'',sampleID) %>% gsub("_downsampled",'',.),
         sampleID = gsub("boreogadus_",'',sampleID)) %>%
  select(sampleID)

# call in some metadata
pop_df <- meta_df %>% # remove everything after space
  left_join(bam_df, ., by = "sampleID") 

########### read in PCA #################################
pca.vectors.list <- list()
varPC1 <- list(); varPC2 <- list()

for(i in 1:length(COVFILES)){
  
  # read in the covariance matrix
  cov <- as.matrix(read.table(COVFILES[i]))
  e <- eigen(cov)                            # calculate eigenvector values
  e_vectors <- as.data.frame(e$vectors) %>% mutate(FID = row_number())  # add FID to e_vectors
  e_per <- e$values/sum(e$values)            # percent explained by each component
  
  # determine the variance explained as a percent
  varPC1[[i]] <- (e$values[1]/sum(e$values))*100 #Variance explained by PC1
  varPC2[[i]] <- (e$values[2]/sum(e$values))*100 #Variance explained by PC2
  
  ##combine row names (population info) with the covariance matrix
  #pca.vectors.list[[i]] = as_tibble(cbind(popFID, e_vectors))[,1:8]
  pca.vectors.list[[i]] = as_tibble(cbind(pop_df, e_vectors))[,1:8]
  
}      

#### PLOT THEME #######################

theme_set(
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(angle = 0, size = 8), axis.title = element_text(size = 10),
        panel.background = element_rect(fill = "white"), panel.spacing = unit(0,"lines"),
        legend.position = "right", strip.text.y = element_text(angle = 0),
        title = element_text(size = 10)
  )
)

mypalette <- color_df$Color2
names(mypalette) <- color_df$Pop

pc12 <- list()

for(j in seq(length(pca.vectors.list))) {  
  
  pca.vectors.i <- as.data.frame(pca.vectors.list[[j]])
  pca.vectors.i$Pop <- factor(pca.vectors.i$Pop, levels = color_df$Pop, ordered = T)
  
  file_name = gsub('boreogadus_','',basename(file_path_sans_ext(COVFILES[j])))
  chr_name = gsub("_.*",'',file_name)
  chr = chrom_df$chr[which(chrom_df$chrName == chr_name)]
  startMb = str_split_i(file_name,"_",2) %>% sub("s","",.) %>% as.numeric(.)/1e6
  endMb = str_split_i(file_name,"_",3) %>% sub("e","",.) %>% as.numeric(.)/1e6
  plot_title = paste('chr',chr,':',startMb,'-',endMb,'Mb')
  
  pc12[[j]] <- ggplot(data = pca.vectors.i, 
                      aes(x=V1, y=V2, fill = Pop)) + 
    geom_point(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
    scale_fill_manual(name = "Location", values = mypalette) +
    ggtitle(plot_title) +
    labs(x = paste0("PC1 (",round(varPC1[[j]], digits = 2),"%)"), 
         y= paste0("PC2 (",round(varPC2[[j]], digits = 2),"%)")) +
    scale_x_continuous(expand = expansion(0.05,0.05)) +
    scale_y_continuous(expand = expansion(0.05,0.05))
  
}

pcs12 <- Reduce('+',pc12) + plot_layout(guides = "collect",
                                        nrow = 2)
pcs12

ggsave(paste0("./figures/pca/20241222/boreogadus_sex",suffix,"_allIndivs.jpg"),
       plot = pcs12, width = 16, height = 8, units = "in")

### Split by Pop ########################################################


pc12 <- list()
pc_facet <- NULL

for(j in seq(length(pca.vectors.list))) { 
  
  pca.vectors.i <- as.data.frame(pca.vectors.list[[j]])
  pca.vectors.i$Pop <- factor(pca.vectors.i$Pop, levels = color_df$Pop, ordered = T)
  
  file_name = gsub('boreogadus_','',basename(file_path_sans_ext(COVFILES[j])))
  
  matches <- str_match(file_name, "_s(\\d+)_e(\\d+)")
  start <- as.numeric(matches[2])/1e6
  end <- as.numeric(matches[3])/1e6
  
  chr = chrom_df$chr[which(chrom_df$chrName == chr_name)]
  
  pc12[[j]] <- ggplot(data = pca.vectors.i, 
                      aes(x=V1, y=V2, fill = Pop)) + 
    geom_point(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
    scale_fill_manual(values = mypalette) +
    ggtitle(paste0('Chr 5: ',start,'-',end,'Mb')) +
    labs(x = paste0("PC1 (",round(varPC1[[j]], digits = 2),"%)"), 
         y= paste0("PC2 (",round(varPC2[[j]], digits = 2),"%)")) +
    scale_x_continuous(expand = expansion(0.05,0.05)) +
    scale_y_continuous(expand = expansion(0.05,0.05))
  
  pc_facet[[j]] <- pca.vectors.i %>%
    filter(V1 > .01 | V1 < -.01) %>%
    ggplot(aes(x=V1, y=V2, fill = Pop)) + 
    geom_jitter(alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
    scale_fill_manual(values = mypalette) +
    facet_wrap(~ Pop,nrow=2) +
    ggtitle(paste0('Chr 5: ',start,'-',end,'Mb')) +
    labs(x = paste0("PC1 (",round(varPC1[[j]], digits = 2),"%)"), 
         y= paste0("PC2 (",round(varPC2[[j]], digits = 2),"%)")) +
    scale_x_continuous(expand = expansion(0.05,0.05), breaks=c(-0.1,0,0.1)) +
    scale_y_continuous(expand = expansion(0.05,0.05)) +
    geom_vline(xintercept=0,linetype="dashed",alpha = 0.7,color="gray40")+
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank())
}

pcsex <- Reduce('|',pc12) + plot_layout(guides = "collect")
pcsex

ggsave(paste0("./figures/pca/20241222/boreogadus_sex_regions",suffix,"_allIndivs.jpg"),
       plot = pcsex, width = 10, height = 5, units = "in")

ggsave(paste0("./figures/pca/20241222/boreogadus_sex_byPop",suffix,"_allIndivs_2026.jpg"),
       plot = pc_facet[[2]], width = 10, height = 5, units = "in")


left <- as.data.frame(pca.vectors.list[[2]]) %>%
  filter(V1 < 0) %>%
  count(Pop) %>% rename(sex1 = n)

right <- as.data.frame(pca.vectors.list[[2]]) %>%
  filter(V1 >= 0) %>%
  count(Pop) %>% rename(sex2 = n)

predicted_sex_count_by_pop <- left_join(left, right)

### Let's incorporate the heterozygosity from the inversion Rmd here ##########

het <- read.table("./results/inversion/het/boreogadus_OZ177908.1_s25700000_e27600000.het", header=T) %>%
        mutate(het = (N.NM. - O.HOM.) / N.NM.,
               rel_het = ( het - min(het,na.rm=T)) / (max(het,na.rm=T) - min(het,na.rm=T)) ) %>%
        bind_cols(bam_df) %>%
  select(sampleID, het, rel_het)

# weird conversions for plotting
pca.variables <- pca.vectors.i %>%
  select(sampleID, Pop, V1, V2) %>%
  merge(., het, by="sampleID") %>%
  mutate(rel_het = as.character(rel_het),
         het = as.character(het)) %>%
  rename(Heterozygosity = het,
         Location = Pop) %>%
  pivot_longer(c(Heterozygosity,Location),
               names_to = "variable",
               values_to = "value")

pca.variables %>%
  filter(V1 > .01 | V1 < -.01) %>%
  ggplot(aes(x=V1, y=V2)) + 
  geom_jitter(aes(fill = value), 
              data = ~ subset(., variable == "Location"),
              alpha = 0.7, size = 2.5, pch = 21, color = 'gray20',
              show.legend = F) +
  scale_fill_manual(name = "Location", values = mypalette) +
  new_scale_fill() +
  geom_jitter(aes(fill=as.numeric(value)),
              data = ~ subset(., variable == "Heterozygosity"),
              alpha = 0.7, size = 2.5, pch = 21, color = 'gray20') +
  scale_fill_viridis_b(name = "Heterozygosity", breaks = seq(0,0.5,by=0.1)) +
  facet_wrap(~ variable,nrow=2) +
  geom_vline(xintercept=0,linetype="dashed",alpha = 0.7,color="gray40")+
  ggtitle(paste0('Chr 5: ',start,'-',end,'Mb')) +
  labs(x = paste0("PC1 (",round(varPC1[[j]], digits = 2),"%)"), 
       y= paste0("PC2 (",round(varPC2[[j]], digits = 2),"%)")) +
  scale_x_continuous(expand = expansion(0.05,0.05), breaks=c(-0.1,0,0.1)) +
  scale_y_continuous(expand = expansion(0.05,0.05)) +
  theme_bw() +
  theme(panel.grid = element_blank())


pc_facet[[2]] + theme(axis.text = element_blank(),
                      axis.title = element_blank())

