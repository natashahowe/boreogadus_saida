# Boreogadus Saida
# Bams lists for FST
# PCA filter lists to calculate PCAs by subsetting
# Natasha Howe

library(tidyverse)
library(stringr)
library(tools)

#### INPUT NAMING AND FILES ####

METADATAFILE <- "./data/R/arctic_cod_pop_metadata.csv"
prefix <- "boreogadus"
bamname <- "filtered" # using bamslist filtered by depth

#### METADATA #####
metadata <- read.csv(METADATAFILE, header = T)

# convert bam to FID with sampleID
bamlist <- read.table(paste0("./data/bam/",prefix,"_",bamname,"_bams.txt"), 
                     header = F) %>%
  dplyr::mutate(bam = V1,
                temp = basename(file_path_sans_ext(V1)),
                temp = gsub("_sorted","",temp),
                sampleID = gsub("_clipped","",temp)) %>%
  select(c(sampleID, bam))

# for if there needs to be removal of individuals for fst
#depths_df <- read.delim("./data/R/boreogadus_depths.csv", sep = "\t", header = F, row.names = NULL, col.names = c("sampleID", "depth"))

df <- left_join(bamlist, metadata2, by = "sampleID")
View(df)

# VIEW SUMMARY OF METADATA 
df %>% group_by(Population) %>% summarize(n=n())


#### BAMS LISTS ######

Labrador <- df %>%
  filter(Population == "Labrador") %>% select(bam)
nrow(Labrador)
write.table(Labrador, "./data/bam/boreogadus_Labrador_bams.txt", 
            col.names=F,row.names=F,quote=F)

Iceland <- df %>%
  filter(Population == "Iceland") %>% select(bam)
nrow(Iceland)
write.table(Iceland, "./data/bam/boreogadus_Iceland_bams.txt", 
            col.names=F,row.names=F,quote=F)

Pond <- df %>%
  filter(Population == "Pond Inlet") %>% select(bam)
nrow(Pond)
write.table(Pond, "./data/bam/boreogadus_Pond_bams.txt", 
            col.names=F,row.names=F,quote=F)

Coronation <- df %>%
  filter(Population == "Coronation Gulf") %>% select(bam)
nrow(Coronation)
write.table(Coronation, "./data/bam/boreogadus_Coronation_bams.txt", 
            col.names=F,row.names=F,quote=F)

Broughton <- df %>%
  filter(Population == "Broughton Island") %>% select(bam)
nrow(Broughton)
write.table(Broughton, "./data/bam/boreogadus_Broughton_bams.txt", 
            col.names=F,row.names=F,quote=F)

Chukchi <- df %>%
  filter(Population == "Chukchi Sea") %>% select(bam)
nrow(Chukchi)
write.table(Chukchi, "./data/bam/boreogadus_Chukchi_bams.txt", 
            col.names=F,row.names=F,quote=F)

Frobisher <- df %>%
  filter(Population == "Frobisher Bay") %>% select(bam)
nrow(Frobisher)
write.table(Frobisher, "./data/bam/boreogadus_Frobisher_bams.txt", 
            col.names=F,row.names=F,quote=F)

SEBaffin <- df %>%
  filter(Population == "Southeast Baffin") %>% select(bam)
nrow(SEBaffin)
write.table(SEBaffin, "./data/bam/boreogadus_SEBaffin_bams.txt", 
            col.names=F,row.names=F,quote=F)

Southhampton <- df %>%
  filter(Population == "Southhampton Is") %>% select(bam)
nrow(Southhampton)
write.table(Southhampton, "./data/bam/boreogadus_Southhampton_bams.txt", 
            col.names=F,row.names=F,quote=F)

Hudson <- df %>%
  filter(Population == "Hudson Strait") %>% select(bam)
nrow(Hudson)
write.table(Hudson, "./data/bam/boreogadus_Hudson_bams.txt", 
            col.names=F,row.names=F,quote=F)

### Subset populations with larger sample sizes

ChukchiSub <- Chukchi[1:18,]
nrow(ChukchiSub)
write.table(ChukchiSub, "./data/bam/boreogadus_Chukchi_subset_bams.txt", 
            col.names=F,row.names=F,quote=F)

PondSub <- Pond[1:18,]
  nrow(PondSub)
  write.table(PondSub,"./data/bam/boreogadus_Pond_subset_bams.txt",
            sep ="\t",col.names=F,row.names=F,quote=F)


  
#### PCA FILTERS ######

nochukPCA <- df %>%
  mutate(Filt = ifelse(Population=="Chukchi Sea",0,1)) %>% select(Filt)
  nochukPCA %>% group_by(Filt) %>% summarise(n())
write.table(nochukPCA, "./data/pca_filter/noChukchi_pcafilter.txt", 
            col.names=F,row.names=F,quote=F)

noChukIcePCA <- df %>%
  mutate(Filt = ifelse(Population=="Chukchi Sea"|Population=="Iceland",0,1)) %>% select(Filt)
noChukIcePCA %>% group_by(Filt) %>% summarise(n())
write.table(noChukIcePCA, "./data/pca_filter/noChukchiIceland_pcafilter.txt", 
            col.names=F,row.names=F,quote=F)

noChukCorPCA <- df %>%
  mutate(Filt = ifelse(Population=="Chukchi Sea"|Population=="Coronation Gulf",0,1)) %>% select(Filt)
noChukCorPCA %>% group_by(Filt) %>% summarise(n())
write.table(noChukCorPCA, "./data/pca_filter/noChukchiCoronation_pcafilter.txt", 
            col.names=F,row.names=F,quote=F)

noIcePCA <- df %>%
  mutate(Filt = ifelse(Population=="Iceland",0,1)) %>% select(Filt)
noIcePCA %>% group_by(Filt) %>% summarise(n())
write.table(noIcePCA, "./data/pca_filter/noIceland_pcafilter.txt", 
            col.names=F,row.names=F,quote=F)
