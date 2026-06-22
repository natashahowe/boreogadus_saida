
# This was example code, but I sent this to Matt for him to plot...
# I think he must have done it a little different

chrom_df <- read.table('./data/R/boreogadus_chromosomes.txt', header = T)


inversions <- read_delim("./data/inversion_sheet.csv") %>%
  mutate(chr = as.numeric(gsub("Chr","",chr.)),
         start.pos = gsub("-.*$","",mb_pos),
         end.pos = sub(".*-","",mb_pos))

inversion <- inversions %>%
  filter(pca_three_clusters == "Yes",
         high_het_prop == "Yes") %>%
  inner_join(chrom_df) %>%
  select(-mb_pos, -chr.) %>%
  select(chr, chrName, start.pos, end.pos, everything())

write.csv(inversion, "./results/inversion/inversion_cleaned_data.csv",
          row.names = F)

test_df <- read.delim2("./results/fst/boreogadus_Chukchi-Pond_20241222.fst.SNP.txt", 
                       header = T, sep = "\t", row.names = NULL, 
                       col.names = c("region", "chrName", "midPos", "Nsites", "Fst")) %>% 
  select(-c(region, Nsites))

# check that both position files are in MB
# I did it by chrName since I tested it with FST, can do it with chr column instead if MDS data is in that format
inversionMDS <- mds_df %>% 
  rowwise() %>%
  mutate(inv_region = any(inversion$chrName == chrName &
                            midPos >= inversion$start.pos &
                            midPos <= inversion$end.pos))

ggplot(inversionMDS) +
  geom_point(aes(x = midPos, y = MDS, color = inv_region)) +
  scale_color_manual(values = c("black","dodgerblue")) +
  theme_minimal() + theme(axis.text.x = element_text(hjust = 0.5)) # hjust will place the chr-text in the middle rather than end
