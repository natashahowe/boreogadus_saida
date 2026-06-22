
packages_needed <- c("tidyverse", "stringr", "gtools", "reshape2", "ggplot2")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i],  repos = "http://cran.us.r-project.org")}
  library(packages_needed[i], character.only = TRUE)
}
library(here)

here()
CHROM="OZ177908.1"
PREFIX <- paste0("boreogadus_",CHROM,"_plink")

### MATT'S LD CODE ###############

matt_ld <- read.delim(paste0("./results/ld/", PREFIX,"_r2.ld"),sep="",header=F,stringsAsFactors = FALSE)
colnames(matt_ld) <- c("BP_A","BP_B","R2")
matt_ld <- matt_ld[order(matt_ld$R2),]
ldplot <- matt_ld %>% 
  filter(R2>=0.3) %>% 
  ggplot(data=.)+
  geom_point(aes(x=BP_A,y=BP_B,color=R2),shape=15,size=0.7)+
  scale_color_gradient()

tiff(paste0("./figures/ld/", PREFIX,"_LD_r2_plink.tiff"), width = 14, height = 12, res = 300, units = "in", compression = "lzw")
ldplot
dev.off()

### Matt's plot code but ngsld #####

r <- read.table(paste0("./results/ld/boreogadus_OZ177904.1_10_r2.ld"), 
                header=FALSE, stringsAsFactors=FALSE)
colnames(r) <- c("snp1", "snp2", "r2") # remove distance from colnames for plink


### By Location ####


# CHROM="OZ177908.1"
# GROUP <- "Chukchi"
# PREFIX <- paste0("boreogadus_",CHROM,"_",GROUP)
# 
# matt_ld <- read.delim(paste0("./results/ld/", PREFIX,"_r2.ld"),sep="",header=F,stringsAsFactors = FALSE)
# colnames(matt_ld) <- c("BP_A","BP_B","R2")
# matt_ld <- matt_ld[order(matt_ld$R2),]
# ldplot <- matt_ld %>% 
#   filter(R2>=0.3) %>% 
#   ggplot(data=.)+
#   geom_point(aes(x=BP_A,y=BP_B,color=R2),shape=15,size=0.7)+
#   scale_color_gradient()
# 
# tiff(paste0("./figures/ld/", PREFIX,"_LD_r2_plink.tiff"), width = 14, height = 12, res = 300, units = "in", compression = "lzw")
# ldplot
# dev.off()
# 



