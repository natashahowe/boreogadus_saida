# ADMIXTURE ANALYSIS

library(tidyverse)
library(reshape2)
library(patchwork)
library(tibble)
library(here) 
library(viridis)
library(stringr)

### Required Input Variables ###################################################
MAX_K = 6 # colors are hard coded for K=6 or less. Need to change if increase!!

# Different runs
runs <- c('wholegenome', 'fused_chroms','fused_adaptive_alpha0.05',
          'unfused_neutral_alpha0.05_e4_rmPCAngsdOutliers_rmInversions_rmchr9a','unfused_adaptive_alpha0.05', 'unfused_chroms')
run_names <- c('Whole Genome', 'Fused: All Loci', 'Fused: Adaptive Loci', 'Unfused: Neutral Loci', 'Unfused: Adaptive Loci', 'Unfused: All Loci')
# NAME='wholegenome'
# NAME='fused_chroms'
# NAME='fused_adaptive_alpha0.05'
# NAME='unfused_neutral_alpha0.05'
# NAME='unfused_adaptive_alpha0.05'

### Data Inputs ################################################################
meta_df <- read.csv( "./data/R/boreogadus_metadata.csv", header = T) %>%
  mutate(Population = ifelse(Population=='Southhampton Island', 'Southampton Island', Population),
         Pop = case_when(Population == 'SE Baffin' ~ 'SEBaffin',
                         TRUE ~ gsub(' .*$','',Population))) %>%
  filter(Pop != 'Unknown')

color_df <- read.delim2(paste0(here(),"./data/R/color_metadata_allpops_flip.txt"), 
                        header = T, row.names = NULL, sep = "\t") %>%
  mutate(Pop = gsub('hh','h',Pop)) %>%
  left_join(meta_df %>% distinct(Population, Pop), by  = "Pop")

BAMFILE <- "./data/R/boreogadus_filtered_downsampled_bams.txt"

bam_df <- read.table(BAMFILE, header = F) %>%
  mutate(sampleID = basename(file_path_sans_ext(V1))%>% 
           gsub("boreogadus_|_downsampled|_sorted", "", .)) %>%
  select(sampleID)

# call in some metadata
pop_df <- meta_df %>% # remove everything after space
  left_join(bam_df, ., by = "sampleID") 

admix_colors <- c("K1" = 'navy',"K2" = 'skyblue2',"K3" = '#2C967D',
                  "K4" = '#FDAA5C',"K5" = 'violetred4',"K6" = 'khaki2')

mypalette <- as.vector(color_df$Color2) # turn colors into vector
  names(mypalette) <- color_df$Pop # attach pop name to palette color
  pop.factor.levels <- color_df$Pop

### RUNS LOOP ##############################################################
for(NAME in runs){
  
### Start loop over K-values ################################################
  admix_plot <- list()
  
  for(K in 2:MAX_K){
      
      WILDCARD = paste0('*_',NAME,'*k',K,"-1*.qopt")
        ADMIXFILES <- Sys.glob(file.path(here(),"results","admix",WILDCARD))
        basename(ADMIXFILES)
    
      admix_ks <- read.delim2(ADMIXFILES[1], sep = "", row.names = NULL, header = F,
                              col.names = paste0('K',1:K))
      
      admix_df <- as.data.frame(lapply(admix_ks, as.numeric))
    
      #### Prepare data for plotting ############
      
      # combine FID and admix output file to denote indiv and pop
      admix_df = as_tibble(cbind(pop_df, admix_df))
      admix_df$Pop <- factor(admix_df$Pop, levels = unique(color_df$Pop))
    
      # make column names say K instead of V
      admix_col_names <- c("Indiv", "Pop_Old", "Population", "Pop", "K1")
      # account for number of Ks in column naming
      for(i in 2:K){
        admix_col_names <- append(admix_col_names, paste0("K",i))
      }
      colnames(admix_df) <- admix_col_names
      
      # change structure of admix_df to plot
      admixed <- admix_df %>%
        select(-Pop_Old) %>%
        melt(., id = c("Indiv", "Population", "Pop"))
      colnames(admixed) <- c("Indiv", "Population", "Pop", "Cluster", "Proportion")
      
      # arrange individuals in population by their max cluster for better viewing
      order_by_cluster <- admixed %>%
        group_by(Indiv, Pop) %>%
        summarise(
          dom_Cluster = Cluster[which.max(Proportion)],
          max_prop = max(Proportion),
          .groups = "drop"
        ) %>%
        group_by(Pop) %>%
        arrange(Pop, dom_Cluster, desc(max_prop)) %>%
        mutate(within_pop_order = row_number()) %>%
        ungroup()
    
      admixed <- admixed %>%
        left_join(order_by_cluster %>% select(Indiv, within_pop_order), by = "Indiv")
      
      pop_levels <- unique(admixed$Pop[order(admixed$Pop)])
      
      admixed <- admixed %>%
        arrange(factor(Pop, levels = pop_levels), within_pop_order) %>%
        mutate(Indiv = factor(Indiv, levels = unique(Indiv)))
      
      #### Plot #################################################
      # population bar
      
      ## Make the colored bar that is for your population colors 
      ## add a column called Locality_text that just has Locality as all the values 
      ## this will serve as your X-axis and your population ID will serve as your y-axis
      admixed <- admixed %>% 
        mutate(Locality_text = "Locality")
      
      # admix data
      admix_plot[[K]] <- ggplot(admixed, aes(x = Indiv, y = Proportion, fill = Cluster)) +
        geom_bar(width=1, stat = "identity", position = "stack") +
        facet_grid( ~ Pop, scales = "free_x", space = "free_x") +
        labs(y = paste0("K",K)) +
        theme_minimal() +
        theme(axis.line = element_blank(), 
              title = element_blank(),
              legend.text = element_text(size = 22),
              legend.title = element_text(size = 25),
              plot.margin = unit(c(t = 0,r = 0,l = 0,b = 0), "cm"),
              axis.text.y = element_text(size = 18, color = 'black'),
              axis.text.x = element_blank(),
              axis.ticks.y = element_line(),
              axis.title.y = element_text(hjust=1, angle=0, size = 25),
              axis.title.x = element_blank(), 
              strip.text.x = element_blank()) +
        scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1),
                           expand = expansion(mult = c(0, 0))) +
        scale_fill_manual(values = admix_colors)
      
    }
    
    # plot a bar that represents population ids for the individuals in rows of the genotype heatmap
    pop_plot <- ggplot(admixed, aes(x = Indiv, y = Locality_text, fill = Pop)) + 
      geom_tile() + 
      facet_grid( ~ Pop, scales = "free_x", space = "free_x", switch = "x") +
      theme_minimal() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(), 
            legend.position = "none", plot.margin = unit(c(0,0,0,0), "cm"),
            axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            title = element_blank(),
            strip.text.x = element_text(size = 22, vjust = 0.5, hjust = 1, angle = 90),
            strip.placement = "inside", strip.background.x = element_blank()) +
      xlab(label = "Population") +
      scale_fill_manual(values = mypalette)  
    
  admix_plot[-length(admix_plot)] <- lapply(admix_plot[-length(admix_plot)], 
                                            function(p) p + theme(legend.position = "none"))
  
  admix_all_Ks <- Reduce('/',Filter(Negate(is.null), admix_plot)) / pop_plot + 
    plot_layout(guides = "collect",
                heights = c(rep(1,length(Filter(Negate(is.null), admix_plot))),0.25)) +
    plot_annotation(title = run_names[which(runs == NAME)],
                    theme = theme(plot.title = element_text(hjust=0.5, size=25,
                                                            margin = margin(0,0,10,0))))
  admix_all_Ks 
  
  ### Max Likelihood K Calculation ###############################################
  # from Bay Lab marine omics webpage
  # need 3 log files per K
  NAME
  
  WILDCARD = paste0('*_',NAME,'*.log')
  LOGFILES <- Sys.glob(file.path(here(),"results","admix",WILDCARD))
  basename(LOGFILES)
  
  nIter = length(LOGFILES[grep("k2",LOGFILES)])
  nIter
  
  log_list <-lapply(1:length(LOGFILES), 
                    FUN = function(i) readLines(LOGFILES[i]))
  
  foundset<-sapply(1:length(LOGFILES), 
                   FUN= function(x) log_list[[x]][which(str_sub(log_list[[x]], 1, 1) == 'b')])
  
  #now lets store it in a dataframe
  #make a dataframe with an index 2:MAX K values, this corresponds to our K values
  logs <-data.frame(K = rep(2:MAX_K, each=nIter),
                    like = as.vector(gsub('best like=','',foundset))) %>% # likelihood values
    mutate(like = gsub("^(-?[0-9]+\\.[0-9]+).*", "\\1",like),
           like = as.numeric(like))
  # head(logs, n = 11)
  # is.numeric(logs$like)
  
  #and now we can calculate our delta K and probability
  deltaK <- tapply(logs$like, logs$K, 
                   FUN= function(x) mean(abs(x))/sd(abs(x)))
  
  maxlikeK <- names(which.max(deltaK))
  
  # deltaK
  # maxlikeK
  # plot with max likelihood mentioned
  jpeg(paste0("./figures/admix/combine/boreogadus_",NAME,"_maxLikeK-",maxlikeK,"_",format(Sys.Date(), "%Y%m%d"),".jpeg"), 
       width = 14, height = 16, res = 150, units = "in")
  print(admix_all_Ks)
  dev.off()
  
  pdf(paste0("./figures/admix/combine/boreogadus_",NAME,"_maxLikeK-",maxlikeK,"_",format(Sys.Date(), "%Y%m%d"),".pdf"), 
       width = 14, height = 16)
  print(admix_all_Ks)
  dev.off()

}
  
  