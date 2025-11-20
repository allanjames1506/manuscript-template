# a meta analysis of Arabidopsis clock ChIP data and overlap with cold TF network

# 1 LIBRARIES----

#library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(janitor)
library(purrr)
library(ggh4x)
library(MetaCycle)
devtools::install_github("vqv/ggbiplot")
library(ggbiplot)
library(ggrepel)
library(circlize)
library(UpSetR)
library(devtools)
#install_github("jokergoo/ComplexHeatmap")
#library(ComplexHeatmap)
library(ComplexUpset)
library(igraph)
library(ndtv)
library(animation)
library(ggbreak)
library(patchwork)
library(readxl)
library(readr)
library(gt)
library(gtExtras)
library(GOplot)
library(org.At.tair.db)
library(clusterProfiler)
library(cowplot)
library(forcats)
library(ggthemes)
library(ggtext)

# 2 CLUSTER GROUP PROFILES - WebPlotDigitizer----
# *2.1 DE profiles----
# gather and join together all the WebPlotDigitizer files for each cluster group
# WebPlotDigitizer https://apps.automeris.io/wpd/
clusters_aggregated_DE <- list.files(path = './01_tidy_data/cluster_image_analysis_aggregate_DE', 
                                  pattern = '*.csv', full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  purrr::reduce(full_join, by = 'id') %>% 
  janitor::remove_empty(which = 'cols') %>% 
  relocate(c(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5, cluster_6, cluster_7, 
             cluster_8, cluster_9), .before = cluster_10) 

# split by individual days
clusters_aggregated_DE_day1 <- clusters_aggregated_DE[row.names(clusters_aggregated_DE) %in% 1:17, ]
clusters_aggregated_DE_day2 <- clusters_aggregated_DE[row.names(clusters_aggregated_DE) %in% 17:33, ]
clusters_aggregated_DE_day5 <- clusters_aggregated_DE[row.names(clusters_aggregated_DE) %in% 34:50, ]

# transpose so that columns are sample numbers (s1 to s9 in half steps) and rows are cluster groups, and rename columns
clusters_aggregated_DE_day1_df <- data.frame(t(clusters_aggregated_DE_day1[, -1])) %>% 
  rename_with(~ paste0("s", seq(from = 1, to = 9, by = 0.5))) %>% 
  mutate(clusters = paste0("cluster_", 1:12)) %>% 
  relocate(clusters) 

clusters_aggregated_DE_day2_df <- data.frame(t(clusters_aggregated_DE_day2[, -1])) %>% 
  rename_with(~ paste0("s", seq(from = 9, to = 17, by = 0.5))) %>% 
  mutate(clusters = paste0("cluster_", 1:12)) %>% 
  relocate(clusters) 

clusters_aggregated_DE_day5_df <- data.frame(t(clusters_aggregated_DE_day5[, -1])) %>% 
  rename_with(~ paste0("s", seq(from = 18, to = 26, by = 0.5))) %>% 
  mutate(clusters = paste0("cluster_", 1:12)) %>% 
  relocate(clusters) 

# **2.1.1 Selected clusters plotted ----

clusters_aggregated_DE_pivot_longer <- clusters_aggregated_DE %>% 
  pivot_longer(cols = starts_with ('cluster'),
               names_to = "cluster", 
               values_to = "z_score") %>%
  mutate(id = case_when(id == 1 ~ 0,
                        id == 1.5 ~ 1.5,
                        id == 2 ~ 3,
                        id == 2.5 ~ 4.5,
                        id == 3 ~ 6,
                        id == 3.5 ~ 7.5,
                        id == 4 ~ 9,
                        id == 4.5 ~ 10.5,
                        id == 5 ~ 12,
                        id == 5.5 ~ 13.5,
                        id == 6 ~ 15,
                        id == 6.5 ~ 16.5,
                        id == 7 ~ 18,
                        id == 7.5 ~ 19.5,
                        id == 8 ~ 21,
                        id == 8.5 ~ 22.5,
                        id == 9 ~ 24,
                        id == 9.5 ~ 25.5,
                        id == 10 ~ 27,
                        id == 10.5 ~ 28.5,
                        id == 11 ~ 30,
                        id == 11.5 ~ 31.5,
                        id == 12 ~ 33,
                        id == 12.5 ~ 34.5,
                        id == 13 ~ 36,
                        id == 13.5 ~ 37.5,
                        id == 14 ~ 39,
                        id == 14.5 ~ 40.5,
                        id == 15 ~ 42,
                        id == 15.5 ~ 43.5,
                        id == 16 ~ 45,
                        id == 16.5 ~ 46.5,
                        id == 17 ~ 48,
                        id == 18 ~ 96,
                        id == 18.5 ~ 97.5,
                        id == 19 ~ 99,
                        id == 19.5 ~ 100.5,
                        id == 20 ~ 102,
                        id == 20.5 ~ 103.5,
                        id == 21 ~ 105,
                        id == 21.5 ~ 106.5,
                        id == 22 ~ 108,
                        id == 22.5 ~ 109.5,
                        id == 23 ~ 111,
                        id == 23.5 ~ 112.5,
                        id == 24 ~ 114,
                        id == 24.5 ~ 115.5,
                        id == 25 ~ 117,
                        id == 25.5 ~ 118.5,
                        id == 26 ~ 120),
         days = case_when(id <= 48 ~ 'Day 1 - Day 2',
                          TRUE ~ 'Day 5'),
         cluster = case_when(cluster == 'cluster_1' ~ 'cluster 1',
                             cluster == 'cluster_2' ~ 'cluster 2',
                             cluster == 'cluster_3' ~ 'cluster 3',
                             cluster == 'cluster_4' ~ 'cluster 4',
                             cluster == 'cluster_5' ~ 'cluster 5',
                             cluster == 'cluster_6' ~ 'cluster 6',
                             cluster == 'cluster_7' ~ 'cluster 7',
                             cluster == 'cluster_8' ~ 'cluster 8',
                             cluster == 'cluster_9' ~ 'cluster 9',
                             cluster == 'cluster_10' ~ 'cluster 10',
                             cluster == 'cluster_11' ~ 'cluster 11',
                             cluster == 'cluster_12' ~ 'cluster 12',
                             TRUE ~ NA))

clusters_aggregated_DE_pivot_longer <- clusters_aggregated_DE_pivot_longer %>% 
  mutate(days = factor(days, levels = c('Day 1 - Day 2', 'Day 5')))

annotations_DE <- data.frame(
  label = c('Day 1', 'Day 2', 'Day 5'),
  cluster = c('cluster 1', 'cluster 1', 'cluster 1'),
  x     = c(12, 44, 108),
  y     = c(2, 2, 2)
)

clusters_plot_DE <- clusters_aggregated_DE_pivot_longer %>%
  ggplot(aes(x=id, y=z_score)) +
  geom_vline(xintercept = 24, col = "lightblue", size = 2) +
  geom_line() +
  geom_point() +
  theme_linedraw() +
  xlim(-1, 122) +
  scale_x_break(c(48, 96)) +
  scale_x_continuous(breaks = seq(0, 120, 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust=0.5, size = 16, face = "bold"),
        axis.text = element_text(size=12),
        axis.title=element_text(size=16,face="bold")) +
  annotate("rect", xmin = c(0, 24, 96), xmax = c(12, 36, 108), ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "grey50") +
  labs(y = "z-score", x = "hours") +
  facet_grid(factor(cluster, levels = c('cluster 1',
                                        'cluster 2',
                                        'cluster 3',
                                        'cluster 4',
                                        'cluster 5',
                                        'cluster 6',
                                        'cluster 7',
                                        'cluster 8',
                                        'cluster 9',
                                        'cluster 10',
                                        'cluster 11',
                                        'cluster 12')) ~ ., scales = "free_y") + 
  theme(strip.text = element_text(face = "bold", color = "grey20", hjust = 0.5, size = 14),
                                                     strip.background = element_rect(fill = "lightblue", linetype = "solid",
                                                                                     color = "black", linewidth = 1)) +
  geom_text(data = annotations_DE,
    mapping = aes(x = x, y = y, label = label),
    size = 6,
    fontface = 'bold'
  )

clusters_plot_DE <- clusters_plot_DE + 
  labs(title = 'DE gene clusters') +
  theme(plot.title = element_text(size=16, color='black', face='bold', hjust =0))

ggsave(clusters_plot_DE, file = './03_plots/clusters_plot_DE1.png', width=6, height=16, units="in",dpi=300)

# *2.2 DTU profiles----
# gather and join together all the WebPlotDigitizer files for each cluster group
# WebPlotDigitizer https://apps.automeris.io/wpd/
clusters_aggregated_DTU <- list.files(path = './01_tidy_data/cluster_image_analysis_aggregate_DTU', 
                                  pattern = '*.csv', full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  purrr::reduce(full_join, by = 'id') %>% 
  janitor::remove_empty(which = 'cols') %>% 
  relocate(c(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5, cluster_6, cluster_7, 
             cluster_8, cluster_9), .before = cluster_10) 

# split by individual days
clusters_aggregated_DTU_day1 <- clusters_aggregated_DTU[row.names(clusters_aggregated_DTU) %in% 1:17, ]
clusters_aggregated_DTU_day2 <- clusters_aggregated_DTU[row.names(clusters_aggregated_DTU) %in% 17:33, ]
clusters_aggregated_DTU_day5 <- clusters_aggregated_DTU[row.names(clusters_aggregated_DTU) %in% 34:50, ]

# transpose so that columns are sample numbers (s1 to s9 in half steps) and rows are cluster groups, and rename columns
clusters_aggregated_DTU_day1_df <- data.frame(t(clusters_aggregated_DTU_day1[, -1])) %>% 
  rename_with(~ paste0("s", seq(from = 1, to = 9, by = 0.5))) %>% 
  mutate(clusters = paste0("cluster_", 1:10)) %>% 
  relocate(clusters) 

clusters_aggregated_DTU_day2_df <- data.frame(t(clusters_aggregated_DTU_day2[, -1])) %>%   rename_with(~ paste0("s", seq(from = 9, to = 17, by = 0.5))) %>% 
  mutate(clusters = paste0("cluster_", 1:10)) %>% 
  relocate(clusters) 

clusters_aggregated_DTU_day5_df <- data.frame(t(clusters_aggregated_DTU_day5[, -1])) %>%   rename_with(~ paste0("s", seq(from = 18, to = 26, by = 0.5))) %>% 
  mutate(clusters = paste0("cluster_", 1:10)) %>% 
  relocate(clusters) 

# **2.2.1 Selected clusters plotted ----

clusters_aggregated_DTU_pivot_longer <- clusters_aggregated_DTU %>% 
  pivot_longer(cols = starts_with ('cluster'),
               names_to = "cluster", 
               values_to = "z_score") %>%
  mutate(id = case_when(id == 1 ~ 0,
                        id == 1.5 ~ 1.5,
                        id == 2 ~ 3,
                        id == 2.5 ~ 4.5,
                        id == 3 ~ 6,
                        id == 3.5 ~ 7.5,
                        id == 4 ~ 9,
                        id == 4.5 ~ 10.5,
                        id == 5 ~ 12,
                        id == 5.5 ~ 13.5,
                        id == 6 ~ 15,
                        id == 6.5 ~ 16.5,
                        id == 7 ~ 18,
                        id == 7.5 ~ 19.5,
                        id == 8 ~ 21,
                        id == 8.5 ~ 22.5,
                        id == 9 ~ 24,
                        id == 9.5 ~ 25.5,
                        id == 10 ~ 27,
                        id == 10.5 ~ 28.5,
                        id == 11 ~ 30,
                        id == 11.5 ~ 31.5,
                        id == 12 ~ 33,
                        id == 12.5 ~ 34.5,
                        id == 13 ~ 36,
                        id == 13.5 ~ 37.5,
                        id == 14 ~ 39,
                        id == 14.5 ~ 40.5,
                        id == 15 ~ 42,
                        id == 15.5 ~ 43.5,
                        id == 16 ~ 45,
                        id == 16.5 ~ 46.5,
                        id == 17 ~ 48,
                        id == 18 ~ 96,
                        id == 18.5 ~ 97.5,
                        id == 19 ~ 99,
                        id == 19.5 ~ 100.5,
                        id == 20 ~ 102,
                        id == 20.5 ~ 103.5,
                        id == 21 ~ 105,
                        id == 21.5 ~ 106.5,
                        id == 22 ~ 108,
                        id == 22.5 ~ 109.5,
                        id == 23 ~ 111,
                        id == 23.5 ~ 112.5,
                        id == 24 ~ 114,
                        id == 24.5 ~ 115.5,
                        id == 25 ~ 117,
                        id == 25.5 ~ 118.5,
                        id == 26 ~ 120),
         days = case_when(id <= 48 ~ 'Day 1 - Day 2',
                          TRUE ~ 'Day 5'),
         cluster = case_when(cluster == 'cluster_1' ~ 'cluster 1',
                             cluster == 'cluster_2' ~ 'cluster 2',
                             cluster == 'cluster_3' ~ 'cluster 3',
                             cluster == 'cluster_4' ~ 'cluster 4',
                             cluster == 'cluster_5' ~ 'cluster 5',
                             cluster == 'cluster_6' ~ 'cluster 6',
                             cluster == 'cluster_7' ~ 'cluster 7',
                             cluster == 'cluster_8' ~ 'cluster 8',
                             cluster == 'cluster_9' ~ 'cluster 9',
                             cluster == 'cluster_10' ~ 'cluster 10',
                             TRUE ~ NA))

clusters_aggregated_DTU_pivot_longer <- clusters_aggregated_DTU_pivot_longer %>% 
  mutate(days = factor(days, levels = c('Day 1 - Day 2', 'Day 5')))

annotations_DTU <- data.frame(
  label = c('Day 1', 'Day 2', 'Day 5'),
  cluster = c('cluster 1', 'cluster 1', 'cluster 1'),
  x     = c(12, 36, 108),
  y     = c(1.2, 1.2, 1.2)
)

clusters_plot_DTU <- clusters_aggregated_DTU_pivot_longer %>%
  ggplot(aes(x=id, y=z_score)) +
  geom_vline(xintercept = 24, col = "lightblue", size = 2) +
  geom_line() +
  geom_point() +
  theme_linedraw() +
  xlim(-1, 122) +
  scale_x_break(c(48, 96)) +
  scale_x_continuous(breaks = seq(0, 120, 6)) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust=0.5, size = 16, face = "bold"),
        axis.text = element_text(size=12),
        axis.title=element_text(size=16,face="bold")) +
  annotate("rect", xmin = c(0, 24, 96), xmax = c(12, 36, 108), ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "grey50") +
  labs(y = "", x = "hours") +
  facet_grid(factor(cluster, levels = c('cluster 1',
                                        'cluster 2',
                                        'cluster 3',
                                        'cluster 4',
                                        'cluster 5',
                                        'cluster 6',
                                        'cluster 7',
                                        'cluster 8',
                                        'cluster 9',
                                        'cluster 10')) ~ ., scales = "free_y") + 
  theme(strip.text = element_text(face = "bold", color = "grey20", hjust = 0.5, size = 14),
        strip.background = element_rect(fill = "lightblue", linetype = "solid",
                                        color = "black", linewidth = 1)) +
  geom_text(data = annotations_DTU,
            mapping = aes(x = x, y = y, label = label),
            size = 6,
            fontface = 'bold'
  )

clusters_plot_DTU <- clusters_plot_DTU + 
  labs(title = 'DTU transcript clusters') +
  theme(plot.title = element_text(size=16, color='black', face='bold', hjust =0))

ggsave(clusters_plot_DTU, file = './03_plots/clusters_plot_DTU1.png', width=6, height=16, units="in",dpi=300)

#*2.3 Plots combined----

combined_DE_DTU <-clusters_plot_DE + (clusters_plot_DTU/plot_spacer() + plot_layout(heights = c(5.5,1))) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = 'bold'))

ggsave(combined_DE_DTU, file = './03_plots/combined_DE_DTU1.png', width=12, height=16, units="in",dpi=300)

ggsave(combined_DE_DTU, file = './03_plots/combined_DE_DTU1.png', width=12, height=16, units="in",dpi=300)


# 1 TF CLUSTERS and GENE_IDs----

# An analysis of established and published Arabidopsis clock ChIP targets in TF network
# read the TF network clusters
# The DE network (Suppl dataset 7A) is 6810 genes over 12 clusters (clusters 1-12)

calixto_S7A <- read_csv("./00_raw_data/Calixto_suppl_dataset_7A_DE.csv") %>%
  select(2:13) %>% 
  pivot_longer(cols = starts_with ('cluster'),
               names_to = "cluster", 
               values_to = "gene_ID",
               values_drop_na = TRUE) %>% 
  filter(grepl('AT', gene_ID)) %>%
  mutate(cluster = str_sub(cluster, 9, -1))

calixto_S7B <- read_csv("./00_raw_data/Calixto_suppl_dataset_7B_DTU.csv") %>%
  select(2:11) %>% 
  pivot_longer(cols = starts_with ('cluster'),
               names_to = "cluster", 
               values_to = "gene_ID",
               values_drop_na = TRUE) %>% 
  filter(grepl('AT', gene_ID)) %>%
  mutate(cluster = str_sub(cluster, 9, -1),
         gene_ID = str_sub(gene_ID, 1, 9))

# 2 CLOCK ChIP DATASETS----

# *2.1 LHY----

# LHY dataset
# read in the Adams LHY paper dataset and skip first 2 lines
# The LHY paper is Adams et al. (2018) New Phytologist 220(3); 897
# supplemental data set (Table S2) 
adams <- read_csv("./00_raw_data/nph15415-sup-0002-tables2.csv",skip=2) %>% 
  dplyr::select(gene_ID = 1) %>%
  filter(!is.na(gene_ID)) %>%
  distinct(gene_ID) %>%
  mutate(clock = 'LHY') 

# *2.2 CCA1----
# **2.2.1 CCA1 Nagel----

# CCA1 Nagel et al. dataset
# read in the Nagel CCA1 paper dataset and skip first 2 lines
# The CCA1 paper is Nagel et al. (2015) PNAS 112(34); E4802
# supplemental data set (Table S1) 
nagel <- read_csv("./00_raw_data/pnas.1513609112.sd01.csv",skip=2) %>% 
  dplyr::select(gene_ID = 10) %>% 
  mutate(gene_ID = str_sub(gene_ID, end = 9)) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'CCA1-Nagel') 

# **2.2.2 CCA1 Kamioka----

# CCA1 Kamioka et al. dataset
# read in the Kamioka CCA1 paper dataset and skip first 2 lines
# The CCA1 paper is Kamioka et al. (2016) Plant Cell 28(3); 696
# supplemental data set (Table S1C)
kamioka <- read_csv("./00_raw_data/TPC2015-00737-RAR3_Supplemental_Data_Set_1C.csv",skip=3) %>% 
  dplyr::select(gene_ID = 10) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'CCA1-Kamioka') 

# **2.2.3 Nagel-Kamioka merge----

# merge the nagel and kamioka CCA1 datasets
# use inner_join from dplyr
# 249 obs.
kamioka_nagel_merge <- inner_join(nagel, kamioka, by = "gene_ID") 

kamioka_nagel_merge <- kamioka_nagel_merge %>% 
  mutate(clock = 'CCA1-Nagel-Kamioka')  

# *2.3 TOC1----

# TOC1 dataset
# read in the Huang TOC1 paper dataset
# The TOC1 paper is Huang et al. (2012) Science 336:75
# supplemental data set (Table S1)
huang <- read_csv("./00_raw_data/Huang TOC1 CHiP TableS1.csv") %>% 
  dplyr::select(gene_ID = 14) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'TOC1') 

# *2.4 PRR5----

# PRR5 dataset
# read in the Nakamichi PRR5 paper dataset
# The PRR5 paper is Nakamichi et al. (2012) PNAS 109:17123
# supplemental data set (Table S3) 
nakamichi <- read_csv("./00_raw_data/Dataset S3 Nakamichi et al PRR5 binding targets PNAS 2012.csv", skip=2) %>% 
  dplyr::select(gene_ID = 3) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'PRR5') 

# *2.5 PRR7----

# PRR7 dataset
# read in the Liu PRR7 paper dataset
# The PRR7 paper is Liu et al. (2013) The Plant Journal 76:101
# supplemental data set (Table S1)
liu <- read_csv("./00_raw_data/Dataset S1 Liu et al PRR7 edit.csv") %>% 
  dplyr::select(gene_ID = 17) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'PRR7') 

# *2.6 LUX----

# LUX dataset
# read in the Ezer EC paper for the LUX dataset (LUX_17 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (LUX_17 tab of Table S6) 
ezer_LUX <- read_csv("./00_raw_data/Ezer et al nplants Suppl Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) %>% 
  mutate(clock = 'LUX') 

# *2.7 ELF3----

# ELF3 dataset
# read in the Ezer EC paper for the ELF3 dataset (ELF3_22 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (ELF3_22 tab of Table S6) 
ezer_ELF3 <- read_csv("./00_raw_data/ELF3_22 Ezer Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'ELF3') 

# *2.8 ELF4----

# ELF4 dataset
# read in the Ezer EC paper for the ELF4 dataset (ELF4_22 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (ELF4_22 tab of Table S6) 
ezer_ELF4 <- read_csv("./00_raw_data/ELF4_22 Ezer Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'ELF4') 

CHIPs <- bind_rows(adams,
                   nagel,
                   kamioka,
                   #kamioka_nagel_merge,
                   huang,
                   nakamichi,
                   liu,
                   ezer_LUX,
                   ezer_ELF3,
                   ezer_ELF4)

CHIPs$clock <- factor(CHIPs$clock, levels = c("CCA1-Nagel", "CCA1-Kamioka", "LHY", "TOC1", "PRR5", "PRR7", "LUX", "ELF3", "ELF4"))
                   
summarise_CHIPs <- CHIPs %>% 
  group_by(clock) %>% 
  summarise(CHIPs_total = n()) %>% 
  ungroup()

# 8 MERGE TF with CLOCK TARGETS----
merge_TF_clock <- function(network, df, clock_id){
  
  merge <- inner_join(network, df, by = 'gene_ID') %>%
    arrange(nchar(cluster), cluster) %>% 
    mutate(clock = {{clock_id}})
  
  return(merge)
}

# DE merge
DE_adams_merge <- merge_TF_clock(calixto_S7A, adams, 'LHY')
DE_nagel_merge <- merge_TF_clock(calixto_S7A, nagel, 'CCA1-Nagel')
DE_kamioka_merge <- merge_TF_clock(calixto_S7A, kamioka, 'CCA1-Kamioka')
DE_kamioka_nagel_merge <- merge_TF_clock(calixto_S7A, kamioka_nagel_merge, 'CCA1-Nagel-Kamioka')
DE_huang_merge <- merge_TF_clock(calixto_S7A, huang, 'TOC1')
DE_nakamichi_merge <- merge_TF_clock(calixto_S7A, nakamichi, 'PRR5')
DE_liu_merge <- merge_TF_clock(calixto_S7A, liu, 'PRR7')
DE_ezer_LUX_merge <- merge_TF_clock(calixto_S7A, ezer_LUX, 'LUX')
DE_ezer_ELF3_merge <- merge_TF_clock(calixto_S7A, ezer_ELF3, 'ELF3')
DE_ezer_ELF4_merge <- merge_TF_clock(calixto_S7A, ezer_ELF4, 'ELF4')

DE_clock <- bind_rows(DE_adams_merge,
                      DE_nagel_merge,
                      DE_kamioka_merge,
                      #DE_kamioka_nagel_merge,
                      DE_huang_merge,
                      DE_nakamichi_merge,
                      DE_liu_merge,
                      DE_ezer_LUX_merge,
                      DE_ezer_ELF3_merge,
                      DE_ezer_ELF4_merge) %>% 
  mutate(type = 'DE')

# DTU merge
DTU_adams_merge <- merge_TF_clock(calixto_S7B, adams, 'LHY')
DTU_nagel_merge <- merge_TF_clock(calixto_S7B, nagel, 'CCA1-Nagel')
DTU_kamioka_merge <- merge_TF_clock(calixto_S7B, kamioka, 'CCA1-Kamioka')
DTU_kamioka_nagel_merge <- merge_TF_clock(calixto_S7B, kamioka_nagel_merge, 'CCA1-Nagel-Kamioka')
DTU_huang_merge <- merge_TF_clock(calixto_S7B, huang, 'TOC1')
DTU_nakamichi_merge <- merge_TF_clock(calixto_S7B, nakamichi, 'PRR5')
DTU_liu_merge <- merge_TF_clock(calixto_S7B, liu, 'PRR7')
DTU_ezer_LUX_merge <- merge_TF_clock(calixto_S7B, ezer_LUX, 'LUX')
DTU_ezer_ELF3_merge <- merge_TF_clock(calixto_S7B, ezer_ELF3, 'ELF3')
DTU_ezer_ELF4_merge <- merge_TF_clock(calixto_S7B, ezer_ELF4, 'ELF4')

DTU_clock <- bind_rows(DTU_adams_merge,
                       DTU_nagel_merge,
                       DTU_kamioka_merge,
                       #DTU_kamioka_nagel_merge,
                       DTU_huang_merge,
                       DTU_nakamichi_merge,
                       DTU_liu_merge,
                       DTU_ezer_LUX_merge,
                       DTU_ezer_ELF3_merge,
                       DTU_ezer_ELF4_merge) %>% 
  mutate(type = 'DTU')

DE_DTU_clock <- bind_rows(DE_clock,
                          DTU_clock)

DE_DTU_clock$clock <- factor(DE_DTU_clock$clock, levels = c("CCA1-Nagel", "CCA1-Kamioka", "LHY", "TOC1", "PRR5", "PRR7", "LUX", "ELF3", "ELF4"))

summarise_merge <- DE_DTU_clock %>% 
  group_by(type, clock) %>% 
  summarise(merge_total = n()) %>% 
  ungroup() %>% 
  left_join(summarise_CHIPs, by = 'clock') %>% 
  mutate(proportion = merge_total/CHIPs_total)

plot_order <- 
  summarise_merge %>% 
  filter(type == "DE") %>% 
  arrange(desc(proportion)) %>% 
  mutate(labels = factor(clock))

summarise_merge <- summarise_merge %>%
  mutate(labels = factor(clock, levels = plot_order$labels, ordered = TRUE)) %>% 
  ## add percentage label with `scales::percent()`
  mutate(perc = scales::percent(merge_total / CHIPs_total, accuracy = .1, trim = FALSE))

col_vec <- c(
  'DE' = '#ef8a62',
  'DTU' = '#999999'
)

ggtext_subtitle <- glue::glue(
  "Gene -
  **<span style='color:{col_vec[[1]]}'>{names(col_vec)[[1]]}</span>** - or transcript -
  **<span style='color:{col_vec[[2]]}'>{names(col_vec)[[2]]}</span>** - clusters"
)

CHIP_merge_plot <- ggplot(summarise_merge, aes(x = labels, y=proportion, fill=type), group = labels) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label = perc), vjust = 1.5,
            position = position_dodge(.9), size = 3, colour = case_when('DE' %in% type ~ 'black', TRUE ~ 'floralwhite')) +
  theme_tufte() +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(angle = 30, size = 14, margin = margin(t = 5, unit = 'mm')),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        plot.subtitle = element_textbox_simple()) +
  #scale_fill_manual(values = c("#ef8a62", "#999999")) +
  scale_fill_manual(values = col_vec) +
  scale_y_continuous(labels = scales::percent)+
  labs(title = 'Commonality of cluster components with clock ChIP targets',
       subtitle = ggtext_subtitle,
       x = '',
       y = '')

summarise_merge_table_DE <- summarise_merge %>% 
  filter(type %in% 'DE') %>% 
  select(transfactor = clock, 
         'cis_targets' = CHIPs_total, 
         'DE_clusters' = merge_total) %>%
  mutate(Epitope = case_when(transfactor == 'CCA1-Nagel' ~ 'GFP-CCA1',
                             transfactor == 'CCA1-Kamioka' ~ 'CCA1-FLAG',
                             transfactor == 'LHY' ~ 'native LHY',
                             transfactor == 'TOC1' ~ 'TOC1-YFP',
                             transfactor == 'PRR5' ~ 'PRR5-FLAG',
                             transfactor == 'PRR7' ~ 'HA-PRR7',
                             transfactor == 'LUX' ~ 'LUX-GFP',
                             transfactor == 'ELF3' ~ 'ELF3-MYC',
                             transfactor == 'ELF4' ~ 'ELF4-HA',
                             TRUE ~ NA),
         transfactor = case_when(transfactor == 'CCA1-Nagel' ~ 'CIRCADIAN CLOCK ASSOCIATED-1 (CCA1)',
                                 transfactor == 'CCA1-Kamioka' ~ 'CIRCADIAN CLOCK ASSOCIATED-1 (CCA1)',
                                 transfactor == 'LHY' ~ 'LATE ELONGATED HYPOCOTYL (LHY)',
                                 transfactor == 'TOC1' ~ 'TIMING OF CAB2 EXPRESSION1 (TOC1)',
                                 transfactor == 'PRR5' ~ 'PSEUDO-RESPONSE REGULATOR5 (PRR5)',
                                 transfactor == 'PRR7' ~ 'PSEUDO-RESPONSE REGULATOR7 (PRR7)',
                                 transfactor == 'LUX' ~ 'LUX ARRYTHMO (LUX)',
                                 transfactor == 'ELF3' ~ 'EARLY FLOWERING3 (ELF3)',
                                 transfactor == 'ELF4' ~ 'EARLY FLOWERING4 (ELF4)',
                                 TRUE ~ NA)) %>% 
  select(1, 4, 2, 3)

summarise_merge_table_DTU <- summarise_merge %>% 
  filter(type %in% 'DTU') %>% 
  select(transfactor = clock, 
         'cis_targets' = CHIPs_total, 
         'DTU_clusters' = merge_total) %>%
  mutate(Epitope = case_when(transfactor == 'CCA1-Nagel' ~ 'GFP-CCA1',
                             transfactor == 'CCA1-Kamioka' ~ 'CCA1-FLAG',
                             transfactor == 'LHY' ~ 'native LHY',
                             transfactor == 'TOC1' ~ 'TOC1-YFP',
                             transfactor == 'PRR5' ~ 'PRR5-FLAG',
                             transfactor == 'PRR7' ~ 'HA-PRR7',
                             transfactor == 'LUX' ~ 'LUX-GFP',
                             transfactor == 'ELF3' ~ 'ELF3-MYC',
                             transfactor == 'ELF4' ~ 'ELF4-HA',
                             TRUE ~ NA),
         transfactor = case_when(transfactor == 'CCA1-Nagel' ~ 'CIRCADIAN CLOCK ASSOCIATED-1 (CCA1)',
                                 transfactor == 'CCA1-Kamioka' ~ 'CIRCADIAN CLOCK ASSOCIATED-1 (CCA1)',
                                 transfactor == 'LHY' ~ 'LATE ELONGATED HYPOCOTYL (LHY)',
                                 transfactor == 'TOC1' ~ 'TIMING OF CAB2 EXPRESSION1 (TOC1)',
                                 transfactor == 'PRR5' ~ 'PSEUDO-RESPONSE REGULATOR5 (PRR5)',
                                 transfactor == 'PRR7' ~ 'PSEUDO-RESPONSE REGULATOR7 (PRR7)',
                                 transfactor == 'LUX' ~ 'LUX ARRYTHMO (LUX)',
                                 transfactor == 'ELF3' ~ 'EARLY FLOWERING3 (ELF3)',
                                 transfactor == 'ELF4' ~ 'EARLY FLOWERING4 (ELF4)',
                                 TRUE ~ NA)) %>% 
  select(1, 4, 2, 3)

summarise_table_final <- summarise_merge_table_DE %>% 
  left_join(summarise_merge_table_DTU, by = c('transfactor', 'Epitope', 'cis_targets')) 
  
summarise_table_final <- summarise_table_final %>% 
  mutate(percent_DE = (DE_clusters / cis_targets) *100,
         percent_DTU = (DTU_clusters / cis_targets) * 100) 

# https://mdigi.tools/lighten-color/#999999

gt_tbl <- gt(summarise_table_final) |>
  tab_header(
    title = md("**Clock ChIP datasets and their overlap with gene and transcript co-expressed clusters**"),
    subtitle = " ") |>
  tab_spanner(
    label = "Common targets",
    columns = c('DE_clusters', 'DTU_clusters')
  ) |>
  tab_spanner(
    label = "ChIP datasets",
    columns = c('Epitope', 'cis_targets')
  ) |>
  tab_footnote(
    footnote = md("Table S1 Nagel *et al.* (2015)"),
    locations = cells_body(columns = 'cis_targets', rows = 1)
  ) |>
  tab_footnote(
    footnote = md("Table S1C Kamioka *et al.* (2016)"),
    locations = cells_body(columns = 'cis_targets', rows = 2)
  ) |>
  tab_footnote(
    footnote = md("Table S2 Adams *et al.* (2018)"),
    locations = cells_body(columns = 'cis_targets', rows = 3)
  ) |>
  tab_footnote(
    footnote = md("Table S1 Huang *et al.* (2012)"),
    locations = cells_body(columns = 'cis_targets', rows = 4)
  ) |>
  tab_footnote(
    footnote = md("Table S3 Nakamichi *et al.* (2012)"),
    locations = cells_body(columns = 'cis_targets', rows = 5)
  ) |>
  tab_footnote(
    footnote = md("Dataset S1 Liu *et al.* (2013)"),
    locations = cells_body(columns = 'cis_targets', rows = 6)
  ) |>
  tab_footnote(
    footnote = md("LUX_17 tab of Table S6 Ezer *et al.* (2017)"),
    locations = cells_body(columns = 'cis_targets', rows = 7)
  ) |>
  tab_footnote(
    footnote = md("ELF3_22 tab of Table S6 Ezer *et al.* (2017)"),
    locations = cells_body(columns = 'cis_targets', rows = 8)
  ) |>
  tab_footnote(
    footnote = md("ELF4_22 tab of Table S6 Ezer *et al.* (2013)"),
    locations = cells_body(columns = 'cis_targets', rows = 9)
  ) |>
  fmt_number(
    columns = 'cis_targets',
    sep_mark = ',',
    decimals = 0
  ) |>
  cols_align(
    align = "center",
    columns = c('cis_targets', 'DE_clusters', 'DTU_clusters')
  ) |>
  cols_label(transfactor = md("*trans* factor"),
             cis_targets = md("*cis* targets"),
             DE_clusters = "DE clusters",
             DTU_clusters = "DTU clusters",
             percent_DE = "% DE",
             percent_DTU = "% DTU") %>% 
  gt_plt_bar_pct(column = percent_DE, scaled = TRUE,
                 labels = TRUE,
                 fill = "#ef8a62", background = "#feb797") %>% 
  gt_plt_bar_pct(column = percent_DTU, scaled = TRUE,
                 labels = TRUE,
                 fill = "#999999", background = "#d6d6d6") |>
  cols_align(
    align = "center",
    columns = c('percent_DE', 'percent_DTU')) |>
  tab_spanner(
    label = "Overlap",
    columns = c('percent_DE', 'percent_DTU'))

gt_tbl
  
  
