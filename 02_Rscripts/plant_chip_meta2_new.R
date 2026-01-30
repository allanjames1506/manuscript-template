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
  theme(text = element_text(family = 'Helvetica'),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust=0.5, size = 18, face = "bold"),
        axis.text = element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        plot.tag = element_text(size = 20, face = 'bold')) +
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
  theme(strip.text = element_text(face = "bold", color = "grey20", hjust = 0.5, size = 16),
                                                     strip.background = element_rect(fill = "lightblue", linetype = "solid",
                                                                                     color = "black", linewidth = 1)) +
  geom_text(data = annotations_DE,
    mapping = aes(x = x, y = y, label = label),
    size = 6,
    fontface = 'bold'
  )

clusters_plot_DE <- clusters_plot_DE + 
  labs(title = 'DE gene clusters',
       tag = 'A') +
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
  theme(text = element_text(family = 'Helvetica'),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust=0.5, size = 18, face = "bold"),
        axis.text = element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        plot.tag = element_text(size = 20, face = 'bold')) +
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
  theme(strip.text = element_text(face = "bold", color = "grey20", hjust = 0.5, size = 16),
        strip.background = element_rect(fill = "lightblue", linetype = "solid",
                                        color = "black", linewidth = 1)) +
  geom_text(data = annotations_DTU,
            mapping = aes(x = x, y = y, label = label),
            size = 6,
            fontface = 'bold'
  )

clusters_plot_DTU <- clusters_plot_DTU + 
  labs(title = 'DTU transcript clusters',
       tag = 'B') +
  theme(plot.title = element_text(size=16, color='black', face='bold', hjust =0))

ggsave(clusters_plot_DTU, file = './03_plots/clusters_plot_DTU1.png', width=6, height=16, units="in",dpi=300)

#*2.3 Plots combined----

combined_DE_DTU <-clusters_plot_DE + (clusters_plot_DTU/plot_spacer() + plot_layout(heights = c(5.5,1))) +
  #plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = 'bold'))

ggsave(combined_DE_DTU, file = './03_plots/combined_DE_DTU1.png', width=12, height=16, units="in",dpi=300)

ggsave(combined_DE_DTU, file = './03_plots/combined_DE_DTU1.png', width=12, height=16, units="in",dpi=300)


# # 1 TF CLUSTERS and GENE_IDs----
# 
# # An analysis of established and published Arabidopsis clock ChIP targets in TF network
# # read the TF network clusters
# # The DE network (Suppl dataset 7A) is 6810 genes over 12 clusters (clusters 1-12)
# 
# calixto_S7A <- read_csv("./00_raw_data/Calixto_suppl_dataset_7A_DE.csv") %>%
#   select(2:13) %>% 
#   pivot_longer(cols = starts_with ('cluster'),
#                names_to = "cluster", 
#                values_to = "gene_ID",
#                values_drop_na = TRUE) %>% 
#   filter(grepl('AT', gene_ID)) %>%
#   mutate(cluster = str_sub(cluster, 9, -1))
# 
# calixto_S7B <- read_csv("./00_raw_data/Calixto_suppl_dataset_7B_DTU.csv") %>%
#   select(2:11) %>% 
#   pivot_longer(cols = starts_with ('cluster'),
#                names_to = "cluster", 
#                values_to = "gene_ID",
#                values_drop_na = TRUE) %>% 
#   filter(grepl('AT', gene_ID)) %>%
#   mutate(cluster = str_sub(cluster, 9, -1),
#          gene_ID = str_sub(gene_ID, 1, 9))
# 
# calixto_S7B_summary <- calixto_S7B %>% 
#   group_by(gene_ID) %>% 
#   summarise(n=n()) %>% 
#   ungroup() %>%
#   group_by(n) %>% 
#   summarise(transcripts = n()) %>% 
#   ungroup()
# 
# calixto_S7B_summary_plot <- ggplot(calixto_S7B_summary, aes(x = n, y=transcripts)) +
#   geom_bar(stat='identity', position='dodge') +
#   geom_text(aes(label = transcripts), hjust = c(1.5, 1.5, 1.5, -0.5, -0.5),
#             size = 4, colour = c('floralwhite', 
#                                  'floralwhite', 
#                                  'floralwhite', 
#                                  'gray20', 
#                                  'gray20')) +
#   theme_tufte() +
#   theme(text = element_text(family = 'Helvetica'),
#         axis.text.x = element_text(size = 14, margin = margin(t = 5, unit = 'mm')),
#         axis.text.y = element_text(size = 14),
#         legend.position = "none",
#         plot.tag = element_text(size = 20, face = 'bold')) +
#   coord_flip() +
#   scale_x_reverse() +
#   labs(title = 'DTU clusters: Number of transcripts per loci',
#        subtitle = '3,972 transcripts from 2,241 loci',
#        x = '',
#        y = 'Loci',
#        tag = 'C')
# 
# calixto_S7B_summary_plot
# 
# combined_DE_DTU_summary <- clusters_plot_DE + (clusters_plot_DTU /calixto_S7B_summary_plot + plot_layout(heights = c(5.5,1)))
# 
# combined_DE_DTU_summary
#  
# ggsave(combined_DE_DTU_summary, file = './03_plots/combined_DE_DTU_summary1.png', width=12, height=16, units="in",dpi=300)

# 3 CLOCK ChIP DATASETS----

# LHY
# LHY dataset
# read in the Adams LHY paper dataset and skip first 2 lines
# The LHY paper is Adams et al. (2018) New Phytologist 220(3); 897
# supplemental data set (Table S2) 
adams <- read_csv("./00_raw_data/nph15415-sup-0002-tables2.csv",skip=2) %>% 
  dplyr::select(gene_ID = 1) %>%
  filter(!is.na(gene_ID)) %>%
  distinct(gene_ID) %>%
  mutate(clock = 'LHY') 

# CCA1
# CCA1 Nagel

# CCA1 Nagel et al. dataset
# read in the Nagel CCA1 paper dataset and skip first 2 lines
# The CCA1 paper is Nagel et al. (2015) PNAS 112(34); E4802
# supplemental data set (Table S1) 
nagel <- read_csv("./00_raw_data/pnas.1513609112.sd01.csv",skip=2) %>% 
  dplyr::select(gene_ID = 10) %>% 
  mutate(gene_ID = str_sub(gene_ID, end = 9)) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'CCA1-Nagel') 

# CCA1 Kamioka

# CCA1 Kamioka et al. dataset
# read in the Kamioka CCA1 paper dataset and skip first 2 lines
# The CCA1 paper is Kamioka et al. (2016) Plant Cell 28(3); 696
# supplemental data set (Table S1C)
kamioka <- read_csv("./00_raw_data/TPC2015-00737-RAR3_Supplemental_Data_Set_1C.csv",skip=3) %>% 
  dplyr::select(gene_ID = 10) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'CCA1-Kamioka') 

# Nagel-Kamioka merge

# merge the nagel and kamioka CCA1 datasets
# use inner_join from dplyr
# 249 obs.
kamioka_nagel_merge <- inner_join(nagel, kamioka, by = "gene_ID") 

kamioka_nagel_merge <- kamioka_nagel_merge %>% 
  mutate(clock = 'CCA1-Nagel-Kamioka')  

# TOC1

# TOC1 dataset
# read in the Huang TOC1 paper dataset
# The TOC1 paper is Huang et al. (2012) Science 336:75
# supplemental data set (Table S1)
huang <- read_csv("./00_raw_data/Huang TOC1 CHiP TableS1.csv") %>% 
  dplyr::select(gene_ID = 14) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'TOC1') 

# PRR5

# PRR5 dataset
# read in the Nakamichi PRR5 paper dataset
# The PRR5 paper is Nakamichi et al. (2012) PNAS 109:17123
# supplemental data set (Table S3) 
nakamichi <- read_csv("./00_raw_data/Dataset S3 Nakamichi et al PRR5 binding targets PNAS 2012.csv", skip=2) %>% 
  dplyr::select(gene_ID = 3) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'PRR5') 

# PRR7

# PRR7 dataset
# read in the Liu PRR7 paper dataset
# The PRR7 paper is Liu et al. (2013) The Plant Journal 76:101
# supplemental data set (Table S1)
liu <- read_csv("./00_raw_data/Dataset S1 Liu et al PRR7 edit.csv") %>% 
  dplyr::select(gene_ID = 17) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'PRR7') 

# LUX

# LUX dataset
# read in the Ezer EC paper for the LUX dataset (LUX_17 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (LUX_17 tab of Table S6) 
ezer_LUX <- read_csv("./00_raw_data/Ezer et al nplants Suppl Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) %>% 
  mutate(clock = 'LUX') 

# ELF3

# ELF3 dataset
# read in the Ezer EC paper for the ELF3 dataset (ELF3_22 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (ELF3_22 tab of Table S6) 
ezer_ELF3 <- read_csv("./00_raw_data/ELF3_22 Ezer Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) %>%
  mutate(clock = 'ELF3') 

# ELF4

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

# MERGE TF with CLOCK TARGETS
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
         transfactor = case_when(transfactor == 'CCA1-Nagel' ~ 'CCA1 (Nagel)',
                                 transfactor == 'CCA1-Kamioka' ~ 'CCA1 (Kamioka)',
                                 TRUE ~ transfactor)) %>%  
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
         transfactor = case_when(transfactor == 'CCA1-Nagel' ~ 'CCA1 (Nagel)',
                                 transfactor == 'CCA1-Kamioka' ~ 'CCA1 (Kamioka)',
                                 TRUE ~ transfactor)) %>% 
  select(1, 4, 2, 3)

summarise_table_final <- summarise_merge_table_DE %>% 
  left_join(summarise_merge_table_DTU, by = c('transfactor', 'Epitope', 'cis_targets')) 

summarise_table_final <- summarise_table_final %>% 
  mutate(percent_DE = (DE_clusters / cis_targets) *100,
         percent_DTU = (DTU_clusters / cis_targets) * 100)

# 4 gtTable----
# https://mdigi.tools/lighten-color/#999999

summarise_table_final |>
  gt() |>
  tab_header(
    title = md("**Clock ChIP targets and their overlap with DE and DTU cluster constituents**"),
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
    footnote = md("CIRCADIAN CLOCK ASSOCIATED-1, Table S1 Nagel *et al.* (2015)"),
    locations = cells_body(columns = c('transfactor', 'cis_targets'), rows = 1)
  ) |>
  tab_footnote(
    footnote = md("CIRCADIAN CLOCK ASSOCIATED-1, Table S1C Kamioka *et al.* (2016)"),
    locations = cells_body(columns = c('transfactor', 'cis_targets'), rows = 2)
  ) |>
  tab_footnote(
    footnote = md("LATE ELONGATED HYPOCOTYL, Table S2 Adams *et al.* (2018)"),
    locations = cells_body(columns = c('transfactor', 'cis_targets'), rows = 3)
  ) |>
  tab_footnote(
    footnote = md("TIMING OF CAB2 EXPRESSION1, Table S1 Huang *et al.* (2012)"),
    locations = cells_body(columns = c('transfactor', 'cis_targets'), rows = 4)
  ) |>
  tab_footnote(
    footnote = md("PSEUDO-RESPONSE REGULATOR5, Table S3 Nakamichi *et al.* (2012)"),
    locations = cells_body(columns = c('transfactor', 'cis_targets'), rows = 5)
  ) |>
  tab_footnote(
    footnote = md("PSEUDO-RESPONSE REGULATOR7, Dataset S1 Liu *et al.* (2013)"),
    locations = cells_body(columns = c('transfactor', 'cis_targets'), rows = 6)
  ) |>
  tab_footnote(
    footnote = md("LUX ARRYTHMO, LUX_17 tab of Table S6 Ezer *et al.* (2017)"),
    locations = cells_body(columns = c('transfactor', 'cis_targets'), rows = 7)
  ) |>
  tab_footnote(
    footnote = md("EARLY FLOWERING3, ELF3_22 tab of Table S6 Ezer *et al.* (2017)"),
    locations = cells_body(columns = c('transfactor', 'cis_targets'), rows = 8)
  ) |>
  tab_footnote(
    footnote = md("EARLY FLOWERING4, ELF4_22 tab of Table S6 Ezer *et al.* (2013)"),
    locations = cells_body(columns = c('transfactor', 'cis_targets'), rows = 9)
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

# 5 Clock ChIP target overlap----

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

calixto_S7B_summary <- calixto_S7B %>% 
  group_by(gene_ID) %>% 
  summarise(n=n()) %>% 
  ungroup() %>%
  group_by(n) %>% 
  summarise(transcripts = n()) %>% 
  ungroup()

calixto_S7B_summary_plot <- ggplot(calixto_S7B_summary, aes(x = n, y=transcripts)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label = transcripts), hjust = c(1.5, 1.5, 1.5, -0.5, -0.5),
            size = 3.5, colour = c('floralwhite', 
                                   'floralwhite', 
                                   'floralwhite', 
                                   'gray20', 
                                   'gray20')) +
  theme_tufte() +
  theme(text = element_text(family = 'Helvetica'),
        axis.text.x = element_text(size = 14, margin = margin(t = 5, unit = 'mm')),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        plot.tag = element_text(size = 20, face = 'bold')) +
  coord_flip() +
  scale_x_reverse() +
  labs(title = 'Number of transcripts per loci',
       subtitle = '3,972 transcripts from 2,241 loci',
       x = '',
       y = 'Loci')

plot_order <- 
  summarise_merge %>% 
  filter(type == "DE") %>% 
  arrange(desc(proportion)) %>% 
  mutate(labels = factor(clock))

summarise_merge_plot <- summarise_merge %>%
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

# Bar plot
CHIP_merge_plot <- ggplot(summarise_merge_plot, aes(x = labels, y=proportion, fill=type), group = labels) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label = perc), vjust = 1.5,
            position = position_dodge(.9), size = 2.5, colour = 'floralwhite') +
  theme_tufte() +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(angle = 30, size = 14, margin = margin(t = 5, unit = 'mm')),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        plot.subtitle = element_textbox_simple()) +
  scale_fill_manual(values = col_vec) +
  scale_y_continuous(labels = scales::percent)+
  labs(title = 'Commonality of cluster components with clock ChIP targets',
       subtitle = ggtext_subtitle,
       x = '',
       y = '')

calixto_S7B_summary_plot / CHIP_merge_plot +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 12, face = 'bold')) &
  theme(plot.background = element_rect(fill = "white", colour = "white"))

# Plots combined
combined_merge_table <- CHIP_merge_plot  / wrap_table(gt_tbl, space = 'fixed')
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = 'bold'))

# 6 MetaCycle - RHYTHMIC SIGNALS----
  
# use meta2d from MetaCycle package to detect rhythmic signals from time-series datasets with multiple methods
# https://cran.r-project.org/web/packages/MetaCycle/MetaCycle.pdf
# https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html
# output files are in specified outdir
  
# DE clusters  
meta2d(infile = './01_tidy_data/clusters_aggregated_DE_day1.csv', filestyle = 'csv', outdir = './00_raw_data/cluster_days_output_day1', 
       timepoints = seq(0, 24, by = 1.5))
  
meta2d(infile = './01_tidy_data/clusters_aggregated_DE_day2.csv', filestyle = 'csv', outdir = './00_raw_data/cluster_days_output_day2', 
       timepoints = seq(0, 24, by = 1.5))
  
meta2d(infile = './01_tidy_data/clusters_aggregated_DE_day5.csv', filestyle = 'csv', outdir = './00_raw_data/cluster_days_output_day5', 
       timepoints = seq(0, 24, by = 1.5))

# read back in the meta2d output file with selected columns
DE_day1_output <- read_csv('./00_raw_data/cluster_days_output_day1/meta2d_clusters_aggregated_DE_day1.csv') %>% 
  select(CycID, LS_BH.Q, meta2d_AMP) %>% 
  dplyr::rename(cluster = CycID, d1_sig = LS_BH.Q,  d1_AMP = meta2d_AMP)

DE_day2_output <- read_csv('./00_raw_data/cluster_days_output_day2/meta2d_clusters_aggregated_DE_day2.csv') %>% 
  select(CycID, LS_BH.Q, meta2d_AMP) %>% 
  dplyr::rename(cluster = CycID, d2_sig = LS_BH.Q,  d2_AMP = meta2d_AMP)

DE_day5_output <- read_csv('./00_raw_data/cluster_days_output_day5/meta2d_clusters_aggregated_DE_day5.csv') %>% 
  select(CycID, LS_BH.Q, meta2d_AMP) %>% 
  dplyr::rename(cluster = CycID, d5_sig = LS_BH.Q,  d5_AMP = meta2d_AMP)

# DTU clusters  
meta2d(infile = './01_tidy_data/clusters_aggregated_DTU_day1.csv', filestyle = 'csv', outdir = './00_raw_data/cluster_days_output_day1', 
       timepoints = seq(0, 24, by = 1.5))

meta2d(infile = './01_tidy_data/clusters_aggregated_DTU_day2.csv', filestyle = 'csv', outdir = './00_raw_data/cluster_days_output_day2', 
       timepoints = seq(0, 24, by = 1.5))

meta2d(infile = './01_tidy_data/clusters_aggregated_DTU_day5.csv', filestyle = 'csv', outdir = './00_raw_data/cluster_days_output_day5', 
       timepoints = seq(0, 24, by = 1.5))

# read back in the meta2d output file with selected columns
DTU_day1_output <- read_csv('./00_raw_data/cluster_days_output_day1/meta2d_clusters_aggregated_DTU_day1.csv') %>% 
  select(CycID, LS_BH.Q, meta2d_AMP) %>% 
  dplyr::rename(cluster = CycID, d1_sig = LS_BH.Q,  d1_AMP = meta2d_AMP)

DTU_day2_output <- read_csv('./00_raw_data/cluster_days_output_day2/meta2d_clusters_aggregated_DTU_day2.csv') %>% 
  select(CycID, LS_BH.Q, meta2d_AMP) %>% 
  dplyr::rename(cluster = CycID, d2_sig = LS_BH.Q,  d2_AMP = meta2d_AMP)

DTU_day5_output <- read_csv('./00_raw_data/cluster_days_output_day5/meta2d_clusters_aggregated_DTU_day5.csv') %>% 
  select(CycID, LS_BH.Q, meta2d_AMP) %>% 
  dplyr::rename(cluster = CycID, d5_sig = LS_BH.Q,  d5_AMP = meta2d_AMP)

# *6.1 compare d1d2----
# **6.1.1 DE d1 vs d2----
# compare d1 vs d2 and set rules for classifying amplitude and rhythm changes between days
# _150 means 1.5 fold up or down difference in amplitude
DE_day1_vs_day2 <- full_join(DE_day1_output, DE_day2_output, by = "cluster") %>% 
  mutate(AMP_change = d2_AMP/d1_AMP *100) %>% 
  mutate(sig_flag = case_when(d1_sig > 0.05 & d2_sig > 0.05 ~ 'nr-nr', 
                              d1_sig <= 0.05 & d2_sig > 0.05 ~ 'r-nr',
                              d1_sig > 0.05 & d2_sig <= 0.05 ~ 'nr-r',
                              d1_sig <= 0.05 & d2_sig <= 0.05 ~ 'r-r',
                              TRUE ~ 'rhythmic')) %>% 
  mutate(amp_flag = case_when(AMP_change <= 66.6  & AMP_change > 33.3 ~ 'lose_medium',
                              AMP_change <= 33.3 ~ 'lose_high',
                              AMP_change >= 150 & AMP_change < 300 ~ 'gain_medium',
                              AMP_change >= 300 ~ 'gain_high',
                              TRUE ~ 'other')) %>% 
  mutate(cluster_id = as.character(paste(1:12))) %>% 
  relocate(cluster_id)

# **6.1.2 DTU d1 vs d2----
# compare d1 vs d2 and set rules for classifying amplitude and rhythm changes between days
# _150 means 1.5 fold up or down difference in amplitude
DTU_day1_vs_day2 <- full_join(DTU_day1_output, DTU_day2_output, by = "cluster") %>% 
  mutate(AMP_change = d2_AMP/d1_AMP *100) %>% 
  mutate(sig_flag = case_when(d1_sig > 0.05 & d2_sig > 0.05 ~ 'nr-nr', 
                              d1_sig <= 0.05 & d2_sig > 0.05 ~ 'r-nr',
                              d1_sig > 0.05 & d2_sig <= 0.05 ~ 'nr-r',
                              d1_sig <= 0.05 & d2_sig <= 0.05 ~ 'r-r',
                              TRUE ~ 'rhythmic')) %>% 
  mutate(amp_flag = case_when(AMP_change <= 66.6  & AMP_change > 33.3 ~ 'lose_medium',
                              AMP_change <= 33.3 ~ 'lose_high',
                              AMP_change >= 150 & AMP_change < 300 ~ 'gain_medium',
                              AMP_change >= 300 ~ 'gain_high',
                              TRUE ~ 'other')) %>% 
  mutate(cluster_id = as.character(paste(1:10))) %>% 
  relocate(cluster_id)

# **6.1.3 DE d1 vs d5----
# compare d1 vs d5 and set rules for classifying amplitude and rhythm changes between days
# _150 means 1.5 fold up or down difference in amplitude
DE_day1_vs_day5 <- full_join(DE_day1_output, DE_day5_output, by = "cluster") %>% 
  mutate(AMP_change = d5_AMP/d1_AMP *100) %>% 
  mutate(sig_flag = case_when(d1_sig > 0.05 & d5_sig > 0.05 ~ 'nr-nr', 
                              d1_sig <= 0.05 & d5_sig > 0.05 ~ 'r-nr',
                              d1_sig > 0.05 & d5_sig <= 0.05 ~ 'nr-r',
                              d1_sig <= 0.05 & d5_sig <= 0.05 ~ 'r-r',
                              TRUE ~ 'rhythmic')) %>% 
  mutate(amp_flag = case_when(AMP_change <= 66.6  & AMP_change > 33.3 ~ 'lose_medium',
                              AMP_change <= 33.3 ~ 'lose_high',
                              AMP_change >= 150 & AMP_change < 300 ~ 'gain_medium',
                              AMP_change >= 300 ~ 'gain_high',
                              TRUE ~ 'other')) %>% 
  mutate(cluster_id = as.character(paste(1:12)),
         colour = case_when(amp_flag %in% 'gain_high' ~ '#e41a1c',
                            amp_flag %in% 'gain_medium' ~ '#377eb8',
                            amp_flag %in% 'lose_high' ~ '#4daf4a',
                            amp_flag %in% 'lose_medium' ~ '#984ea3',
                            amp_flag %in% 'gain_high' ~ '#ff7f00',
                            amp_flag %in% 'other' ~ '#bababa',
                            TRUE ~ NA)) %>% 
  relocate(cluster_id)

DE_day1_vs_day5$amp_flag <- factor(DE_day1_vs_day5$amp_flag, levels = c("lose_high", "lose_medium", "other"))

levels(DE_day1_vs_day5$amp_flag)

# **6.1.4 DTU d1 vs d5----
# compare d1 vs d5 and set rules for classifying amplitude and rhythm changes between days
# _150 means 1.5 fold up or down difference in amplitude
DTU_day1_vs_day5 <- full_join(DTU_day1_output, DTU_day5_output, by = "cluster") %>% 
  mutate(AMP_change = d5_AMP/d1_AMP *100) %>% 
  mutate(sig_flag = case_when(d1_sig > 0.05 & d5_sig > 0.05 ~ 'nr-nr', 
                              d1_sig <= 0.05 & d5_sig > 0.05 ~ 'r-nr',
                              d1_sig > 0.05 & d5_sig <= 0.05 ~ 'nr-r',
                              d1_sig <= 0.05 & d5_sig <= 0.05 ~ 'r-r',
                              TRUE ~ 'rhythmic')) %>% 
  mutate(amp_flag = case_when(AMP_change <= 66.6  & AMP_change > 33.3 ~ 'lose_medium',
                              AMP_change <= 33.3 ~ 'lose_high',
                              AMP_change >= 150 & AMP_change < 300 ~ 'gain_medium',
                              AMP_change >= 300 ~ 'gain_high',
                              TRUE ~ 'other')) %>% 
  mutate(cluster_id = as.character(paste(1:10))) %>% 
  relocate(cluster_id)

# 7 MetaCycle SCATTER PLOTS----

# *7.1 compare cluster amplitudes----
# **7.1.1 DE----
# Scatter plot of the MetaCycle outputs for Amplitude coloured by the amp_flag grouping
DE_plot_day1_vs_day2_AMP <- DE_day1_vs_day2 %>% 
  #filter(pval_flag == 'd1_nr_d2_r' | pval_flag == 'd1_r_d2_nr' | pval_flag == 'd1_r_d2_r') %>%
  ggplot(aes(x = d1_AMP, y = d2_AMP, colour = amp_flag, shape = sig_flag)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  scale_x_continuous(limits = c(0, 1.8), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  geom_point(size = 2.5) +
  ggpubr::theme_pubr() +
  theme(legend.position = "bottom", 
        legend.key = element_blank(),
        legend.box.background = element_rect(color = "black"),
        legend.box.margin = margin(t = 1, l = 1),
        plot.title = element_text(color = "grey30", face = 'bold', hjust = 0.5),
        plot.subtitle = element_text(color = "grey30", hjust = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 0.5, face = 'bold'),
        axis.title.x = element_text(face = 'bold')) +
  scale_colour_brewer(palette = "Set1", labels = c("gain - high", "gain - medium", "lose - high", "lose - medium", "other")) +
  geom_text_repel(aes(label = cluster_id), 
                  show.legend = FALSE,
                  max.overlaps = nrow(DE_day1_vs_day2)) +
  labs(colour = "Amplitude", shape = "Rhythm",
       y= "day 2 \namplitude", x = "day 1 amplitude") +
  ggtitle("Day 1 vs Day 2",
          subtitle = "transition from 20C to 4C") +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 0, nrow = 2))

DE_plot_day1_vs_day5_AMP <- DE_day1_vs_day5 %>% 
  #filter(pval_flag == 'd1_nr_d2_r' | pval_flag == 'd1_r_d2_nr' | pval_flag == 'd1_r_d2_r') %>%
  ggplot(aes(x = d1_AMP, y = d5_AMP, colour = amp_flag, shape = sig_flag)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  scale_x_continuous(limits = c(0, 1.8), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  geom_point(size = 2.5) +
  ggpubr::theme_pubr() +
  theme(legend.position = "bottom", 
        legend.key = element_blank(),
        legend.box.background = element_rect(color = "black"),
        legend.box.margin = margin(t = 1, l = 1),
        plot.title = element_text(color = "grey30", face = 'bold', hjust = 0.5),
        plot.subtitle = element_text(color = "grey30", hjust = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 0.5, face = 'bold'),
        axis.title.x = element_text(face = 'bold')) +
  scale_colour_brewer(palette = "Set1", labels = c("lose - high", "lose - medium", "other")) +
  geom_text_repel(aes(label = cluster_id), 
                  show.legend = FALSE,
                  max.overlaps = nrow(DE_day1_vs_day2)) +
  labs(colour = "Amplitude", shape = "Rhythm",
       y= "day 5 \namplitude", x = "day 1 amplitude") +
  ggtitle("Day 1 vs Day 5",
          subtitle = "acclimation to 4C") +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 0, nrow = 2))

# **7.1.2 DTU----
# Scatter plot of the MetaCycle outputs for Amplitude coloured by the amp_flag grouping
DTU_plot_day1_vs_day2_AMP <- DTU_day1_vs_day2 %>% 
  #filter(pval_flag == 'd1_nr_d2_r' | pval_flag == 'd1_r_d2_nr' | pval_flag == 'd1_r_d2_r') %>%
  ggplot(aes(x = d1_AMP, y = d2_AMP, colour = amp_flag, shape = sig_flag)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  scale_x_continuous(limits = c(0, 1.8), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  geom_point(size = 2.5) +
  ggpubr::theme_pubr() +
  theme(legend.position = "bottom", 
        legend.key = element_blank(),
        legend.box.background = element_rect(color = "black"),
        legend.box.margin = margin(t = 1, l = 1),
        plot.title = element_text(color = "grey30", face = 'bold', hjust = 0.5),
        plot.subtitle = element_text(color = "grey30", hjust = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 0.5, face = 'bold'),
        axis.title.x = element_text(face = 'bold')) +
  scale_colour_brewer(palette = "Set1", labels = c("gain - high", "gain - medium", "lose - high", "lose - medium", "other")) +
  geom_text_repel(aes(label = cluster_id), 
                  show.legend = FALSE,
                  max.overlaps = nrow(DE_day1_vs_day2)) +
  labs(colour = "Amplitude", shape = "Rhythm",
       y= "day 2 \namplitude", x = "day 1 amplitude") +
  ggtitle("Day 1 vs Day 2",
          subtitle = "transition from 20C to 4C") +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 0, nrow = 2))

DTU_plot_day1_vs_day5_AMP <- DTU_day1_vs_day5 %>% 
  #filter(pval_flag == 'd1_nr_d2_r' | pval_flag == 'd1_r_d2_nr' | pval_flag == 'd1_r_d2_r') %>%
  ggplot(aes(x = d1_AMP, y = d5_AMP, colour = amp_flag, shape = sig_flag)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  scale_x_continuous(limits = c(0, 1.8), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  geom_point(size = 2.5) +
  ggpubr::theme_pubr() +
  theme(legend.position = "bottom", 
        legend.key = element_blank(),
        legend.box.background = element_rect(color = "black"),
        legend.box.margin = margin(t = 1, l = 1),
        plot.title = element_text(color = "grey30", face = 'bold', hjust = 0.5),
        plot.subtitle = element_text(color = "grey30", hjust = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 0.5, face = 'bold'),
        axis.title.x = element_text(face = 'bold')) +
  scale_colour_brewer(palette = "Set1", labels = c("gain - high", "gain - medium", "lose - high", "lose - medium", "other")) +
  geom_text_repel(aes(label = cluster_id), 
                  show.legend = FALSE,
                  max.overlaps = nrow(DE_day1_vs_day2)) +
  labs(colour = "Amplitude", shape = "Rhythm",
       y= "day 5 \namplitude", x = "day 1 amplitude") +
  ggtitle("Day 1 vs Day 5",
          subtitle = "acclimation to 4C") +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 0, nrow = 2))
  
  
