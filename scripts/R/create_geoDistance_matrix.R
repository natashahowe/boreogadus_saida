# matrix of pairwise fst and matrix of distances
# ISOLATION BY DISTANCE, FIRST GET ALL GEODISTANCES
# FST/(1 - FST)

library(vegan)
library(tidyverse)
library(reshape2)

fst_geo <- read.delim("./results/ibd/fst_geo_2024.txt", sep = "\t", header = T) %>% select(-fst)
ice_geo <- read.csv("./results/ibd/Iceland_geodistance.csv", header = T) 

order_pops <- c("Chukchi","Coronation","Pond","Broughton","SEBaffin","Hudson",
                "Frobisher","Southhampton","Labrador","Iceland")

geototal <- rbind(fst_geo,ice_geo) %>%
  mutate(POP1 = sub("\\..*","",comparison),
         POP2 = sub("^.*\\.","",comparison)) %>% select(POP1,POP2,geo)

geototal$POP1 <- factor(geototal$POP1, levels = order_pops)
geototal$POP2 <- factor(geototal$POP2, levels = order_pops)

geototal_rev <- geototal %>%
  rename(POP1 = POP2, POP2 = POP1) %>% select(POP1,POP2,geo)

geo_matrix <- reshape2::acast(rbind(geototal,geototal_rev), POP1~POP2, value.var = "geo")

geo_matrix <- geo_matrix[row.names(geo_matrix) != 'Labrador',]
geo_matrix <- geo_matrix[,colnames(geo_matrix) != 'Labrador']

write.csv(geo_matrix,"./results/ibd/geodistances_matrix_20250224.csv",row.names = T)

#fst_matrix <- 
 
  
melt(mat)[seq(from=N+1,to=N^2,by=2*(N+1)),]

matrixConvert()
#https://www.rdocumentation.org/packages/otuSummary/versions/0.1.2/topics/matrixConvert



#fst_geo <- fst_geo[-c(1:7),]

# pearsons correlation coeff.
cor(fst_geo$invFST, fst_geo$geo)
round(rval^2, digits = 2)
  
reg <- ggplot(data = fst_geo, aes(x = geo, y = invFST)) +
  geom_point(pch = 21, fill = "dodgerblue", color = "gray30", 
             size = 3, alpha = 0.9) + 
  annotate("text", x = 300, y = 0.041, label = "paste(italic(r), \" = 0.81\")", parse = TRUE) + 
  annotate("text", x = 300, y = 0.043, label = "paste(italic(r) ^ 2, \" = 0.66\")", parse = TRUE) +#label = paste(expression(r2)," = ",round(r2, digits = 2))) + 
  ylab("Fst/(1-Fst)") + xlab("Marine Distance (km)") +
  geom_smooth(method = "lm", se = F, col = "gray30") + 
  theme_minimal()
reg

ggsave(file = "./figures/arctic_cod_ibd_r0.81_20250128.jpeg", 
       plot = reg, width = 8, height = 6)
