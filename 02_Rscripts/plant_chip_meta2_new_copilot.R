
# ============================================================
# 1. LIBRARIES
summarise_table_final |>
  gt() |>
  tab_header(
    title = md("**Clock ChIP targets and their overlap with DE and DTU cluster constituents**"),
    subtitle = " "
  ) |>
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
  cols_label(
    transfactor = md("*trans* factor"),
    cis_targets = md("*cis* targets"),
    DE_clusters = "DE clusters",
    DTU_clusters = "DTU clusters",
    percent_DE = "% DE",
    percent_DTU = "% DTU"
  ) |>
  # gt_plt_bar_pct(
  #   column = percent_DE,
  #   scaled = TRUE,
  #   labels = TRUE,
  #   fill = "#ef8a62",
  #   background = "#feb797"
  # ) |>
  # gt_plt_bar_pct(
  #   column = percent_DTU,
  #   scaled = TRUE,
  #   labels = TRUE,
  #   fill = "#999999",
  #   background = "#d6d6d6"
  # ) |>
  cols_align(
    align = "center",
    columns = c('percent_DE', 'percent_DTU')
  ) |>
  tab_spanner(
    label = "Overlap",
    columns = c('percent_DE', 'percent_DTU')
  )

# ============================================================
# 10. DTU SUMMARY PLOT
# ============================================================

calixto_S7B_summary <- calixto_S7B %>%
  group_by(gene_ID) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(n) %>%
  summarise(transcripts = n()) %>%
  ungroup() %>%
  mutate(
    text_hjust  = if_else(n <= 3, 1.5, -0.5),
    text_colour = if_else(n <= 3, "floralwhite", "gray20")
  )

calixto_S7B_summary_plot <- ggplot(calixto_S7B_summary, aes(x = n, y = transcripts)) +
  geom_col() +
  geom_text(
    aes(label = transcripts, hjust = text_hjust, colour = text_colour),
    size = 3.5
  ) +
  scale_colour_identity() +
  theme_tufte() +
  theme(
    text = element_text(family = 'Helvetica'),
    axis.text.x = element_text(size = 14, margin = margin(t = 5, unit = 'mm')),
    axis.text.y = element_text(size = 14),
    legend.position = "none",
    plot.tag = element_text(size = 20, face = 'bold')
  ) +
  coord_flip() +
  scale_x_reverse() +
  labs(
    title = 'Number of transcripts per loci',
    subtitle = '3,972 transcripts from 2,241 loci',
    x = '',
    y = 'Loci'
  )

# ============================================================
# 11. CLOCK TARGET OVERLAP BAR PLOT
# ============================================================

plot_levels <- summarise_merge %>%
  mutate(clock = case_when(clock == 'CCA1-Kamioka' ~ 'CCA1-K',
                           clock == 'CCA1-Nagel' ~ 'CCA1-N',
                           TRUE ~ clock)) %>% 
  filter(type == "DE") %>%
  arrange(desc(proportion)) %>%
  pull(clock)

summarise_merge_plot <- summarise_merge %>%
  mutate(clock = case_when(clock == 'CCA1-Kamioka' ~ 'CCA1-K',
                           clock == 'CCA1-Nagel' ~ 'CCA1-N',
                           TRUE ~ clock)) %>% 
  mutate(labels = factor(clock, levels = plot_levels),
         perc = scales::percent(proportion, accuracy = 0.1))

col_vec <- c(DE = '#ef8a62', DTU = '#999999')

ggtext_subtitle <- glue::glue(
  "Gene – <b><span style='color:{col_vec['DE']}'>{'DE'}</span></b> or transcript – 
   <b><span style='color:{col_vec['DTU']}'>{'DTU'}</span></b> clusters"
)

# Bar plot
CHIP_merge_plot <- ggplot(summarise_merge_plot, aes(x = labels, y = proportion, fill = type)) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = perc),
    vjust = 1.5,
    position = position_dodge(width = 0.9),
    size = 2.5,
    colour = 'floralwhite'
  ) +
  theme_tufte() +
  theme(
    text = element_text(family = 'Arial'),
    #axis.text.x = element_text(angle = 30, size = 14, margin = margin(t = 5, unit = 'mm')),
    axis.text.x = element_text(angle = 30, size = 14, hjust = 0.5, vjust = 1),
    #axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none",
    plot.subtitle = element_textbox_simple()
  ) +
  scale_fill_manual(values = col_vec) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = 'Commonality of cluster components with clock ChIP targets',
    subtitle = ggtext_subtitle,
    x = '',
    y = ''
  )

# Combined plot
combined_DE_DTU_no2 <- (calixto_S7B_summary_plot / CHIP_merge_plot) +
  plot_annotation(tag_levels = 'A') &
  theme(
    plot.tag = element_text(size = 12, face = 'bold'),
    plot.background = element_rect(fill = "white", colour = "white")
  )

ggsave(combined_DE_DTU_no2, file = './03_plots/copilot/combined_DE_DTU_no2.png', width=7, height=8, units="in",dpi=300)


# ============================================================
# 6. MetaCycle – rhythmic signals (tidy, modular version)
# ============================================================

# Helper: run MetaCycle and read back the summary
run_meta <- function(infile, outdir, prefix) {
  meta2d(
    infile    = infile,
    filestyle = "csv",
    outdir    = outdir,
    timepoints = seq(0, 24, by = 1.5)
  )
  
  read_csv(file.path(outdir, paste0("meta2d_", prefix, ".csv"))) %>%
    select(CycID, LS_BH.Q, meta2d_AMP) %>%
    rename(
      cluster = CycID,
      sig     = LS_BH.Q,
      AMP     = meta2d_AMP
    )
}

# Helper: compare two days for a given type (DE/DTU)
compare_days <- function(d1, d2, n_clusters, type, day_from, day_to) {
  left_join(d1, d2, by = "cluster", suffix = c("_from", "_to")) %>%
    mutate(
      AMP_change = AMP_to / AMP_from * 100,
      sig_flag = case_when(
        sig_from > 0.05 & sig_to > 0.05 ~ "nr-nr",
        sig_from <= 0.05 & sig_to > 0.05 ~ "r-nr",
        sig_from > 0.05 & sig_to <= 0.05 ~ "nr-r",
        sig_from <= 0.05 & sig_to <= 0.05 ~ "r-r",
        TRUE ~ NA
      ),
      amp_flag = case_when(
        AMP_change <= 33.3 ~ "lose high",
        AMP_change <= 66.6 ~ "lose medium",
        AMP_change >= 300  ~ "gain high",
        AMP_change >= 150  ~ "gain medium",
        TRUE ~ "other"
      ),
      cluster_id = as.character(seq_len(n_clusters)),
      type       = type,
      day_from   = day_from,
      day_to     = day_to,
      comparison = paste0(day_from, "_vs_", day_to)
    ) %>%
    mutate(
      amp_flag = factor(
        amp_flag,
        levels = c("gain high", "gain medium", "lose high", "lose medium", "other"),
        labels = c("gain high", "gain medium", "lose high", "lose medium", "other")
      ),
      sig_flag = factor(
        sig_flag,
        levels = c("nr-nr", "r-nr", "nr-r", "r-r")
      )
    )
}

# ------------------------------------------------------------
# 6.1 Run MetaCycle for DE and DTU clusters
# ------------------------------------------------------------

# DE
DE_d1 <- run_meta(
  infile = "./01_tidy_data/clusters_aggregated_DE_day1.csv",
  outdir = "./00_raw_data/cluster_days_output_day1",
  prefix = "clusters_aggregated_DE_day1"
)

DE_d2 <- run_meta(
  infile = "./01_tidy_data/clusters_aggregated_DE_day2.csv",
  outdir = "./00_raw_data/cluster_days_output_day2",
  prefix = "clusters_aggregated_DE_day2"
)

DE_d5 <- run_meta(
  infile = "./01_tidy_data/clusters_aggregated_DE_day5.csv",
  outdir = "./00_raw_data/cluster_days_output_day5",
  prefix = "clusters_aggregated_DE_day5"
)

# DTU
DTU_d1 <- run_meta(
  infile = "./01_tidy_data/clusters_aggregated_DTU_day1.csv",
  outdir = "./00_raw_data/cluster_days_output_day1",
  prefix = "clusters_aggregated_DTU_day1"
)

DTU_d2 <- run_meta(
  infile = "./01_tidy_data/clusters_aggregated_DTU_day2.csv",
  outdir = "./00_raw_data/cluster_days_output_day2",
  prefix = "clusters_aggregated_DTU_day2"
)

DTU_d5 <- run_meta(
  infile = "./01_tidy_data/clusters_aggregated_DTU_day5.csv",
  outdir = "./00_raw_data/cluster_days_output_day5",
  prefix = "clusters_aggregated_DTU_day5"
)

# ------------------------------------------------------------
# 6.2 Combined summary table of all DE/DTU transitions
# ------------------------------------------------------------

DE_day1_vs_day2  <- compare_days(DE_d1,  DE_d2,  n_clusters = 12, type = "DE",  day_from = "d1", day_to = "d2")
DE_day1_vs_day5  <- compare_days(DE_d1,  DE_d5,  n_clusters = 12, type = "DE",  day_from = "d1", day_to = "d5")
DTU_day1_vs_day2 <- compare_days(DTU_d1, DTU_d2, n_clusters = 10, type = "DTU", day_from = "d1", day_to = "d2")
DTU_day1_vs_day5 <- compare_days(DTU_d1, DTU_d5, n_clusters = 10, type = "DTU", day_from = "d1", day_to = "d5")

rhythm_summary_all <- bind_rows(
  DE_day1_vs_day2,
  DE_day1_vs_day5,
  DTU_day1_vs_day2,
  DTU_day1_vs_day5
) %>%
  mutate(
    comparison = recode(comparison,
                        "d1_vs_d2" = "d1 vs d2",
                        "d1_vs_d5" = "d1 vs d5"),
    type       = factor(type, levels = c("DE", "DTU")),
    comparison = factor(comparison, levels = c("d1 vs d2", "d1 vs d5"))
  )

# Optional: write combined summary to disk
#write_csv(rhythm_summary_all, "./01_tidy_data/rhythm_summary_all_DE_DTU.csv")

# ============================================================
# 7. Facet scatterplots of amplitude changes (DE + DTU)
# ============================================================

amp_palette <- c(
  "gain high"   = "#006400",
  "gain medium" = "#66c2a5",
  "lose high"   = "#8b0000",
  "lose medium" = "#fc8d62",
  "other"       = "#bdbdbd"
)

shape_palette <- c(
  "nr-nr"    = 16,  # both arrhythmic
  "r-nr"     = 17,  # lost rhythm
  "nr-r"     = 15,  # gained rhythm
  "r-r"      = 18  # rhythmic both days
)

scatter_amp_facets <- rhythm_summary_all %>%
  ggplot(aes(x = AMP_from, y = AMP_to, colour = amp_flag, shape = sig_flag)) +
  
  # Quadrant shading
  geom_polygon(
    data = data.frame(
      x = c(0, max(rhythm_summary_all$AMP_from), max(rhythm_summary_all$AMP_from)),
      y = c(0, max(rhythm_summary_all$AMP_from), 0)
    ),
    aes(x, y),
    inherit.aes = FALSE,
    fill = "#fee0d2",   # light red (loss zone)
    alpha = 0.3
  ) +
  geom_polygon(
    data = data.frame(
      x = c(0, 0, max(rhythm_summary_all$AMP_from)),
      y = c(0, max(rhythm_summary_all$AMP_to), max(rhythm_summary_all$AMP_to))
    ),
    aes(x, y),
    inherit.aes = FALSE,
    fill = "#e5f5e0",   # light green (gain zone)
    alpha = 0.3
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  
  # Points + labels
  geom_point(size = 4, alpha = 0.75) +
  geom_text_repel(aes(label = cluster_id), show.legend = FALSE) +
  
  facet_grid(type ~ comparison) +
  scale_colour_manual(values = amp_palette, drop = FALSE) +
  scale_shape_manual(values = shape_palette, drop = FALSE) +
  ggpubr::theme_pubr() +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.spacing = unit(10, "pt"),
    legend.box.margin = margin(t = 1, r = 1, b = 1, l = 1),
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.2, "cm"),
    legend.key = element_blank(),
    legend.box.background = element_rect(color = "black"),
    legend.spacing.y = unit(2.5, 'pt'),
    plot.title = element_text(color = "grey30", face = "bold", hjust = 0.5),
    plot.subtitle = element_text(color = "grey30", hjust = 0.5),
    axis.title.y = element_text(angle = 90, vjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
    plot.tag = element_text(size = 18, face = "bold")) +
  labs(
    colour = "Amplitude",
    shape  = "Rhythm",
    x      = "Amplitude (day 1)",
    y      = "Amplitude (day 2 / day 5)",
    title  = "Amplitude and rhythmicity changes across days",
    subtitle = "DE and DTU clusters: d1 vs d2 and d1 vs d5",
    tag = "A"
  ) +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 0, nrow = 2))

# Save panel
ggsave(
  "./03_plots/copilot/rhythm_amp_facets_DE_DTU.png",
  scatter_amp_facets,
  dpi   = 300,
  width = 8,
  height = 10,
  units = "in"
)

# cluster pattern summary plots

summarise_amp <- function(df, label) {
  df %>%
    mutate(amp_simple = case_when(
      grepl("gain", amp_flag) ~ "gain",
      grepl("lose", amp_flag) ~ "lose",
      TRUE ~ "other"
    )) %>%
    group_by(amp_simple) %>%
    summarise(number = n(), .groups = "drop") %>%
    mutate(
      percent = number / sum(number) * 100,
      condition = label
  )
}

amp_df <- bind_rows(
  summarise_amp(DE_day1_vs_day2,  "DE_d1d2"),
  summarise_amp(DE_day1_vs_day5,  "DE_d1d5"),
  summarise_amp(DTU_day1_vs_day2, "DTU_d1d2"),
  summarise_amp(DTU_day1_vs_day5, "DTU_d1d5")
)

summarise_rhythm <- function(df, label) {
  df %>%
    mutate(rhy_simple = case_when(
      sig_flag == "nr-r" ~ "nr-r",
      sig_flag == "r-nr" ~ "r-nr",
      sig_flag == "r-r"  ~ "r-r",
      sig_flag == "nr-nr" ~ "nr-nr",
      TRUE ~ NA
    )) %>%
    group_by(rhy_simple) %>%
    summarise(number = n(), .groups = "drop") %>%
    mutate(
      percent = number / sum(number) * 100,
      condition = label
    )
}

rhy_df <- bind_rows(
  summarise_rhythm(DE_day1_vs_day2,  "DE_d1d2"),
  summarise_rhythm(DE_day1_vs_day5,  "DE_d1d5"),
  summarise_rhythm(DTU_day1_vs_day2, "DTU_d1d2"),
  summarise_rhythm(DTU_day1_vs_day5, "DTU_d1d5")
)

amp_palette_simple <- c(
  gain  = "#006400",
  lose  = "#8b0000",
  other = "#bdbdbd"
)

rhy_palette <- c(
  "nr-r"   = "#66c2a5",
  "r-nr"   = "#fc8d62",
  "r-r"    = "#006400",
  "nr-nr"  = "#bdbdbd"
)

amp_summary_plot <- amp_df %>%
  mutate(
    amp_simple = factor(amp_simple, levels = c("gain", "lose", "other")),
    type = ifelse(grepl("DE", condition), "DE", "DTU"),
    comparison = recode(
      condition,
      "DE_d1d2"  = "d1 vs d2",
      "DE_d1d5"  = "d1 vs d5",
      "DTU_d1d2" = "d1 vs d2",
      "DTU_d1d5" = "d1 vs d5"
    )
  ) %>%
  ggplot(aes(x = comparison, y = percent, fill = amp_simple)) +
  geom_col(width = 0.7) +
  facet_grid(type ~ .) +
  scale_fill_manual(values = amp_palette_simple) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent)),
    position = position_stack(vjust = 0.5),
    size = 4,
    colour = "floralwhite"
  ) +
  labs(fill = "Amplitude", y = "Percent", x = "", tag = "B") +
  theme_pubr() +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.box.spacing = unit(2.5, "pt"),
    legend.box.background = element_rect(color = "black"),
    plot.tag = element_text(size = 18, face = "bold")
  )

ggsave("./03_plots/copilot/amp_summary_plot.png", amp_summary_plot, width = 9.5, height = 5, dpi = 300)


rhy_summary_plot <- rhy_df %>%
  mutate(
    type = ifelse(grepl("DE", condition), "DE", "DTU"),
    comparison = recode(condition,
                        "DE_d1d2" = "d1 vs d2",
                        "DE_d1d5" = "d1 vs d5",
                        "DTU_d1d2" = "d1 vs d2",
                        "DTU_d1d5" = "d1 vs d5"
    )
  ) %>%
  ggplot(aes(x = comparison, y = percent, fill = rhy_simple)) +
  geom_col(width = 0.7) +
  facet_grid(type ~ .) +
  scale_fill_manual(values = rhy_palette) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent)),
    position = position_stack(vjust = 0.5),
    size = 4,
    colour = "floralwhite"
    ) +
  labs(fill = "Rhythm", y = "Percent", x = "", tag = "C") +
  theme_pubr() +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.box.spacing = unit(2.5, "pt"),
    legend.box.background = element_rect(color = "black"),
    plot.tag = element_text(size = 18, face = "bold")
  )

ggsave("./03_plots/copilot/rhy_summary_plot.png", rhy_summary_plot, width = 9.5, height = 5, dpi = 300)

# option B scatter and summary - two columns
# need to take out tag levels in labs and theme

summary_panel <- cowplot::plot_grid(
  amp_summary_plot,
  rhy_summary_plot,
  labels = c("B", "C"),
  ncol = 1
)

final_figure <- cowplot::plot_grid(
  scatter_amp_facets,
  summary_panel,
  labels = c("A", ""),
  ncol = 3,
  rel_widths = c(2, 1)
)

ggsave("./03_plots/copilot/final_amp_rhythm_figure.png", final_figure, width = 19, height = 10, dpi = 300)



# clusters stratified for clock targets----
# Helper: extract clusters for all amplitude categories
get_all_amp_clusters <- function(df) {
  df %>%
    dplyr::transmute(
      cluster      = cluster_id,
      amp_category = as.character(amp_flag)
    )
}

# Helper: Join amplitude clusters to the unified clock table
join_clock_targets <- function(clock_df, amp_df, label_prefix) {
  clock_df %>%
    dplyr::inner_join(amp_df, by = "cluster") %>%
    dplyr::mutate(type = paste0(label_prefix, "_", amp_category))
}

# Helper: bind amplitude categories per clock target
bind_clock_amp_categories <- function(df) {
  df %>%
    dplyr::group_by(clock) %>%
    dplyr::summarise(data = list(dplyr::pick(everything())), .groups = "drop")
}

# Generate amplitude cluster sets for each comparison
DE_d1d2_amp <- get_all_amp_clusters(DE_day1_vs_day2)
DE_d1d5_amp <- get_all_amp_clusters(DE_day1_vs_day5)

DTU_d1d2_amp <- get_all_amp_clusters(DTU_day1_vs_day2)
DTU_d1d5_amp <- get_all_amp_clusters(DTU_day1_vs_day5)

DE_day1_vs_day5 %>% dplyr::count(amp_flag)

# Join amplitude clusters to the unified clock tables
DE_d1d2_clock <- join_clock_targets(DE_clock, DE_d1d2_amp, "DE_d1d2")
DE_d1d5_clock <- join_clock_targets(DE_clock, DE_d1d5_amp, "DE_d1d5")

DTU_d1d2_clock <- join_clock_targets(DTU_clock, DTU_d1d2_amp, "DTU_d1d2")
DTU_d1d5_clock <- join_clock_targets(DTU_clock, DTU_d1d5_amp, "DTU_d1d5")

# Bind amplitude categories per clock target
DE_d1d2_clock_bound <- bind_clock_amp_categories(DE_d1d2_clock)
DE_d1d5_clock_bound <- bind_clock_amp_categories(DE_d1d5_clock)

DTU_d1d2_clock_bound <- bind_clock_amp_categories(DTU_d1d2_clock)
DTU_d1d5_clock_bound <- bind_clock_amp_categories(DTU_d1d5_clock)

get_all_rhy_clusters <- function(df) {
  df %>%
    dplyr::transmute(
      cluster      = cluster_id,
      rhy_category = as.character(sig_flag)
    )
}

# Join rhythmicity clusters to the unified clock table
join_clock_rhythm <- function(clock_df, rhy_df, label_prefix) {
  clock_df %>%
    dplyr::inner_join(rhy_df, by = "cluster") %>%
    dplyr::mutate(type = paste0(label_prefix, "_", rhy_category))
}

# Bind rhythmicity categories per clock target
bind_clock_rhy_categories <- function(df) {
  df %>%
    dplyr::group_by(clock) %>%
    dplyr::summarise(
      data = list(dplyr::pick(everything())),
      .groups = "drop"
    )
}

# DE rhythmicity
DE_d1d2_rhy <- get_all_rhy_clusters(DE_day1_vs_day2)
DE_d1d5_rhy <- get_all_rhy_clusters(DE_day1_vs_day5)

DE_d1d2_clock_rhy <- join_clock_rhythm(DE_clock, DE_d1d2_rhy, "DE_d1d2")
DE_d1d5_clock_rhy <- join_clock_rhythm(DE_clock, DE_d1d5_rhy, "DE_d1d5")

DE_d1d2_clock_rhy_bound <- bind_clock_rhy_categories(DE_d1d2_clock_rhy)
DE_d1d5_clock_rhy_bound <- bind_clock_rhy_categories(DE_d1d5_clock_rhy)

# DTU rhythmicity
DTU_d1d2_rhy <- get_all_rhy_clusters(DTU_day1_vs_day2)
DTU_d1d5_rhy <- get_all_rhy_clusters(DTU_day1_vs_day5)

DTU_d1d2_clock_rhy <- join_clock_rhythm(DTU_clock, DTU_d1d2_rhy, "DTU_d1d2")
DTU_d1d5_clock_rhy <- join_clock_rhythm(DTU_clock, DTU_d1d5_rhy, "DTU_d1d5")

DTU_d1d2_clock_rhy_bound <- bind_clock_rhy_categories(DTU_d1d2_clock_rhy)
DTU_d1d5_clock_rhy_bound <- bind_clock_rhy_categories(DTU_d1d5_clock_rhy)

# Generate all summaries
# Amplitude
# Summarise amplitude per clock target
summarise_clock_amp <- function(df_bound) {
  df_bound %>%
    tidyr::unnest(data) %>%
    dplyr::count(clock, amp_category) %>%
    dplyr::group_by(clock) %>%
    dplyr::mutate(percent = n / sum(n) * 100) %>%
    dplyr::ungroup()
}

DE_d1d2_amp_clock_summary  <- summarise_clock_amp(DE_d1d2_clock_bound)
DE_d1d5_amp_clock_summary  <- summarise_clock_amp(DE_d1d5_clock_bound)
DTU_d1d2_amp_clock_summary <- summarise_clock_amp(DTU_d1d2_clock_bound)
DTU_d1d5_amp_clock_summary <- summarise_clock_amp(DTU_d1d5_clock_bound)

# Rhythmicity
# Summarise rhythmicity per clock target
summarise_clock_rhy <- function(df_bound) {
  df_bound %>%
    tidyr::unnest(data) %>%
    dplyr::count(clock, rhy_category) %>%
    dplyr::group_by(clock) %>%
    dplyr::mutate(percent = n / sum(n) * 100) %>%
    dplyr::ungroup()
}

DE_d1d2_rhy_clock_summary  <- summarise_clock_rhy(DE_d1d2_clock_rhy_bound)
DE_d1d5_rhy_clock_summary  <- summarise_clock_rhy(DE_d1d5_clock_rhy_bound)
DTU_d1d2_rhy_clock_summary <- summarise_clock_rhy(DTU_d1d2_clock_rhy_bound)
DTU_d1d5_rhy_clock_summary <- summarise_clock_rhy(DTU_d1d5_clock_rhy_bound)

# Plot amplitude summary
amp_clock_plot <- function(df, title, show_legend = TRUE) {
  
  amp_levels <- c("gain high", "gain medium", "other", "lose medium", "lose high")
  clock_levels <- c("CCA1-K", "CCA1-N", "LHY", "TOC1", "PRR5", "PRR7", "LUX", "ELF3", "ELF4")
  
  p <- df %>%
    mutate(
      clock = recode(clock,
                     "CCA1-Kamioka" = "CCA1-K",
                     "CCA1-Nagel"   = "CCA1-N"),
      amp_category = factor(amp_category, levels = amp_levels),
      clock        = factor(clock,        levels = clock_levels)
    ) %>%
    ggplot(aes(x = clock, y = percent, fill = amp_category)) +
    geom_col(position = "stack") +
    scale_fill_manual(
      values = c(
        "gain high"   = "#006400",
        "gain medium" = "#66c2a5",
        "other"       = "#bdbdbd",
        "lose medium" = "#fc8d62",
        "lose high"   = "#8b0000"
      ),
      name = "Amplitude change"
    ) +
    labs(title = title, x = "", y = "Percent of targets") +
    theme_pubr() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    )
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

# Plot rhythmicity summary
rhy_clock_plot <- function(df, title, show_legend = TRUE) {
  
  rhy_levels <- c("nr-nr", "r-nr", "nr-r", "r-r")
  clock_levels <- c("CCA1-K", "CCA1-N", "LHY", "TOC1", "PRR5", "PRR7", "LUX", "ELF3", "ELF4")  
  
  p <- df %>%
    mutate(
      clock = recode(clock,
                     "CCA1-Kamioka" = "CCA1-K",
                     "CCA1-Nagel"   = "CCA1-N"),
      rhy_category = factor(rhy_category, levels = rhy_levels),
      clock        = factor(clock,        levels = clock_levels)
    ) %>%
    ggplot(aes(x = clock, y = percent, fill = rhy_category)) +
    geom_col(position = "stack") +
    scale_fill_manual(
      values = c(
        "nr-nr"    = "#bdbdbd",
        "r-nr"     = "#fc8d62",
        "nr-r"     = "#66c2a5",
        "r-r"      = "#006400"
      ),
      name = "Rhythmicity change"
    ) +
    labs(title = title, x = "", y = "Percent of targets") +
    theme_pubr() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    )
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

make_four_panel <- function(p1, p2, p3, p4, 
                            title = "Clock-target behaviour") {
  (p1 | p2) /
    (p3 | p4) +
    patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = 18, face = "bold"))
    )
}

# for the legend
# Build a dummy amplitude dataset
amp_dummy_df <- tibble(
  clock = factor("dummy", levels = "dummy"),
  percent = 1,
  amp_category = factor(
    c("gain high", "gain medium", "other", "lose medium", "lose high"),
    levels = c("gain high", "gain medium", "other", "lose medium", "lose high")
  )
)

# Build a dummy amplitude plot (legend only)
amp_dummy_plot <- amp_clock_plot(amp_dummy_df, "dummy", show_legend = TRUE)

extract_legend <- function(p) {
  cowplot::get_legend(p + theme(legend.position = "right"))
}

amp_legend <- extract_legend(amp_dummy_plot) 

# Build panels WITHOUT legends
pA_amp <- amp_clock_plot(DE_d1d2_amp_clock_summary,  "DE d1→d2", show_legend = FALSE)
pB_amp <- amp_clock_plot(DTU_d1d2_amp_clock_summary, "DTU d1→d2", show_legend = FALSE)
pC_amp <- amp_clock_plot(DE_d1d5_amp_clock_summary,  "DE d1→d5", show_legend = FALSE)
pD_amp <- amp_clock_plot(DTU_d1d5_amp_clock_summary, "DTU d1→d5", show_legend = FALSE)

# Combine the panels and attach the dummy legend
amp_four_panel_shared_legend <- (
  (
    (pA_amp | pB_amp) /
      (pC_amp | pD_amp)
  ) |
    amp_legend
) +
  patchwork::plot_layout(widths = c(4, 1)) +
  patchwork::plot_annotation(
    title = "Amplitude changes across cold transitions",
    tag_levels = "A"
  )

amp_four_panel_shared_legend

ggsave("./03_plots/copilot/amp_four_panel_shared_legend.png", amp_four_panel_shared_legend, width = 8, height = 8, dpi = 300)

# Build a dummy rhythmicity dataset (contains all 5 categories)
rhy_dummy_df <- tibble(
  clock = factor("dummy", levels = "dummy"),
  percent = 1,
  rhy_category = factor(
    c("nr-nr", "r-nr", "nr-r", "r-r"),
    levels = c("nr-nr", "r-nr", "nr-r", "r-r")
  )
)

# Build a dummy rhythmicity plot and extract the legend
rhy_dummy_plot <- rhy_clock_plot(rhy_dummy_df, "dummy", show_legend = TRUE)
rhy_legend <- extract_legend(rhy_dummy_plot)

# Build the four real rhythmicity panels (no legends)
pA_rhy <- rhy_clock_plot(DE_d1d2_rhy_clock_summary,  "DE d1→d2", show_legend = FALSE)
pB_rhy <- rhy_clock_plot(DTU_d1d2_rhy_clock_summary, "DTU d1→d2", show_legend = FALSE)
pC_rhy <- rhy_clock_plot(DE_d1d5_rhy_clock_summary,  "DE d1→d5", show_legend = FALSE)
pD_rhy <- rhy_clock_plot(DTU_d1d5_rhy_clock_summary, "DTU d1→d5", show_legend = FALSE)

#Assemble the 4‑panel rhythmicity figure with shared legend
rhy_four_panel_shared_legend <- (
  (
    (pA_rhy | pB_rhy) /
      (pC_rhy | pD_rhy)
  ) |
    rhy_legend
) +
  patchwork::plot_layout(widths = c(4, 1)) +
  patchwork::plot_annotation(
    title = "Rhythmicity changes across cold transitions",
    tag_levels = "A"
  ) 

rhy_four_panel_shared_legend

ggsave("./03_plots/copilot/rhy_four_panel_shared_legend.png", rhy_four_panel_shared_legend, width = 8, height = 8, dpi = 300)


# Cluster size plot
# Cluster size plot (simple, grey, sorted, labels inside bars)
DE_cluster_sizes <- DE_DTU_clock %>%
  filter(type == 'DE') %>% 
  dplyr::count(cluster) %>%
  mutate(cluster = reorder(cluster, -n),
         vjust_pos = case_when(
           cluster %in% c("6", "1", "9") ~ 1.2,     # inside the bar
           TRUE           ~ -0.3     # above the bar) 
           ))# highest → lowest

DTU_cluster_sizes <- DE_DTU_clock %>%
  filter(type == 'DTU') %>% 
  dplyr::count(cluster) %>%
  mutate(cluster = reorder(cluster, -n),
         vjust_pos = case_when(
           cluster == "1" ~ 1.2,     # inside the bar
           TRUE           ~ -0.3     # above the bar) 
         ))# highest → lowest

DE_cluster_order <- DE_cluster_sizes %>%
  arrange(desc(n)) %>%
  pull(cluster)

DTU_cluster_order <- DTU_cluster_sizes %>%
  arrange(desc(n)) %>%
  pull(cluster)

DE_cluster_size_plot <- ggplot(DE_cluster_sizes, aes(x = cluster, y = n)) +
  geom_col(fill = "grey70") +
  geom_text(aes(label = n, vjust = vjust_pos), size = 3.5) +
  labs(x = "Cluster", y = "Total genes") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5, 5, 5, 5)
  ) + 
  labs(tag = "A")

DTU_cluster_size_plot <- ggplot(DTU_cluster_sizes, aes(x = cluster, y = n)) +
  geom_col(fill = "grey70") +
  geom_text(aes(label = n, vjust = vjust_pos), size = 3.5) +
  labs(x = "Cluster", y = "Total genes") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5, 5, 5, 5)
  ) + 
  labs(tag = "C")

# Phase‑strip annotation
# Circadian phase strip (with labels + border + title)
phase_df <- tibble(
  phase = factor(c("Dawn", "Day", "Dusk", "Night"),
                 levels = c("Dawn", "Day", "Dusk", "Night")),
  x = 1:4
)

phase_palette <- c(
  "Dawn"  = "#238b45",
  "Day"   = "#fdae61",
  "Dusk"  = "#4d4d4d",
  "Night" = "#8073ac"
)

phase_strip <- ggplot(phase_df, aes(x = x, y = 1, fill = phase)) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_text(aes(label = phase), color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = phase_palette) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Circadian phase") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 12, face = "bold"),
    plot.margin = margin(0, 0, 2, 0),
    legend.position = 'none',
    plot.tag = element_blank()
  )

# Linear stacked bar plot clock targets per cluster group----
# Amplitude----
# *DE----
# Define the biological stack order
clock_order <- c(
  "CCA1-Nagel", "CCA1-Kamioka", "LHY",   # dawn / morning
  "PRR5", "PRR7",                         # daytime
  "TOC1",                                 # evening
  "LUX", "ELF3", "ELF4"                   # night
)

# Join two data frames
DE_df_d1d2 <- DE_DTU_clock %>%
  filter(type == 'DE') %>% 
  left_join(DE_d1d2_amp, by = "cluster") %>%
  mutate(dataset = "d1d2")

DE_df_d1d5 <- DE_DTU_clock %>%
  filter(type == 'DE') %>%
  left_join(DE_d1d5_amp, by = "cluster") %>%
  mutate(dataset = "d1d5")

# Count genes per cluster × clock gene
DE_df_d1d2_counts <- DE_df_d1d2 %>%
  dplyr::count(amp_category, cluster, clock)

DE_df_d1d5_counts <- DE_df_d1d5 %>%
  dplyr::count(amp_category, cluster, clock)

# Compute total counts per cluster (for sorting)
DE_cluster_totals_d1d2 <- DE_df_d1d2_counts %>%
  group_by(amp_category, cluster) %>%
  summarise(total = sum(n), .groups = "drop")

DE_cluster_totals_d1d5 <- DE_df_d1d5_counts %>%
  group_by(amp_category, cluster) %>%
  summarise(total = sum(n), .groups = "drop")

# Merge totals back and sort clusters within each facet
DE_df_d1d2_plot <- DE_df_d1d2_counts %>%
  left_join(DE_cluster_totals_d1d2, by = c("amp_category", "cluster")) %>%
  mutate(
    clock = factor(clock, levels = clock_order),
    cluster = factor(cluster, levels = DE_cluster_order),
    dataset = 'DE_d1d2'
  )

DE_df_d1d5_plot <- DE_df_d1d5_counts %>%
  left_join(DE_cluster_totals_d1d5, by = c("amp_category", "cluster")) %>%
  mutate(
    clock = factor(clock, levels = clock_order),
    cluster = factor(cluster, levels = DE_cluster_order),
    dataset = 'DE_d1d5'
  )

DE_df_all <- bind_rows(DE_df_d1d2_plot, DE_df_d1d5_plot) %>% 
  mutate(dataset = case_when(dataset == 'DE_d1d2' ~ 'DE d1d2',
                             dataset == 'DE_d1d5' ~ 'DE d1d5',
                             TRUE ~ dataset))

# Build the faceted stacked bar plot
tf_palette <- c(
  "CCA1-Nagel"   = "#006d2c",  # deep forest green
  "CCA1-Kamioka" = "#238b45",  # mid forest green
  "LHY"          = "#5aae61",  # your existing light green
  
  "PRR5"         = "#fdae61",  # warm daytime orange
  "PRR7"         = "#f46d43",  # stronger late‑day orange
  
  "TOC1"         = "#4d4d4d",  # dusk grey
  
  "LUX"          = "#542788",  # deep night purple
  "ELF3"         = "#8073ac",  # mid purple
  "ELF4"         = "#b2abd2"   # pale night purple
)

DE_stacked_plot <- ggplot(DE_df_all, aes(x = cluster, y = n, fill = clock)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = tf_palette, breaks = clock_order) +
  scale_x_reordered() +
  facet_grid(amp_category ~ dataset, scales = "free_x") +
  labs(
    x = "Cluster",
    y = "Gene count",
    fill = "Clock context"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box.spacing = unit(2.5, "pt"),
    legend.box.background = element_rect(color = "black")
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) + 
  labs(tag = "B")

DE_final_figure <- 
  DE_cluster_size_plot /
  phase_strip /
  DE_stacked_plot +
  plot_layout(heights = c(0.1, 0.03, 0.87)) 

ggsave("./03_plots/copilot/DE_clock_targets_in_clusters.png", DE_final_figure, width = 6, height = 12, dpi = 300)

# *DTU----
# Join two data frames
DTU_df_d1d2 <- DE_DTU_clock %>%
  filter(type == 'DTU') %>% 
  left_join(DTU_d1d2_amp, by = "cluster") %>%
  mutate(dataset = "d1d2")

DTU_df_d1d5 <- DE_DTU_clock %>%
  filter(type == 'DTU') %>%
  left_join(DTU_d1d5_amp, by = "cluster") %>%
  mutate(dataset = "d1d5")

# Count genes per cluster × clock gene
DTU_df_d1d2_counts <- DTU_df_d1d2 %>%
  dplyr::count(amp_category, cluster, clock)

DTU_df_d1d5_counts <- DTU_df_d1d5 %>%
  dplyr::count(amp_category, cluster, clock)

# Compute total counts per cluster (for sorting)
DTU_cluster_totals_d1d2 <- DTU_df_d1d2_counts %>%
  group_by(amp_category, cluster) %>%
  summarise(total = sum(n), .groups = "drop")

DTU_cluster_totals_d1d5 <- DTU_df_d1d5_counts %>%
  group_by(amp_category, cluster) %>%
  summarise(total = sum(n), .groups = "drop")

# Merge totals back and sort clusters within each facet
DTU_df_d1d2_plot <- DTU_df_d1d2_counts %>%
  left_join(DTU_cluster_totals_d1d2, by = c("amp_category", "cluster")) %>%
  mutate(
    clock = factor(clock, levels = clock_order),
    cluster = factor(cluster, levels = DTU_cluster_order),
    dataset = 'DTU_d1d2'
  )

DTU_df_d1d5_plot <- DTU_df_d1d5_counts %>%
  left_join(DTU_cluster_totals_d1d5, by = c("amp_category", "cluster")) %>%
  mutate(
    clock = factor(clock, levels = clock_order),
    cluster = factor(cluster, levels = DTU_cluster_order),
    dataset = 'DTU_d1d5'
  )

DTU_df_all <- bind_rows(DTU_df_d1d2_plot, DTU_df_d1d5_plot) %>% 
  mutate(dataset = case_when(dataset == 'DTU_d1d2' ~ 'DTU d1d2',
                             dataset == 'DTU_d1d5' ~ 'DTU d1d5',
                             TRUE ~ dataset))

DTU_stacked_plot <- ggplot(DTU_df_all, aes(x = cluster, y = n, fill = clock)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = tf_palette, breaks = clock_order) +
  scale_x_reordered() +
  facet_grid(amp_category ~ dataset, scales = "free_x") +
  labs(
    x = "Cluster",
    y = "Gene count",
    fill = "Clock context"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box.spacing = unit(2.5, "pt"),
    legend.box.background = element_rect(color = "black")
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) + 
  labs(tag = "D")

DTU_final_figure <- 
  DTU_cluster_size_plot /
  phase_strip /
  DTU_stacked_plot +
  plot_layout(heights = c(0.1, 0.03, 0.87)) 

ggsave("./03_plots/copilot/DTU_clock_targets_in_clusters.png", DTU_final_figure, width = 6, height = 12, dpi = 300)

# patched plot - amp----
DE_stacked_noleg <- DE_stacked_plot + theme(legend.position = "none")
DTU_stacked_noleg <- DTU_stacked_plot + theme(legend.position = "none")

shared_legend_f <- extract_legend(DE_stacked_plot)

top_row <- 
  (DE_cluster_size_plot + labs(tag = "A")) |
  (DTU_cluster_size_plot + labs(tag = "C"))

bottom_row <- 
  (DE_stacked_noleg + labs(tag = "B")) |
  (DTU_stacked_noleg + labs(tag = "D"))

grid_2x2 <- 
  top_row /
  phase_strip /
  bottom_row +
  plot_layout(widths = c(1,1), heights = c(0.2, 0.03, 1)) 

grid_2x2_legend <- 
  wrap_plots(grid_2x2) /
  plot_spacer() /
  shared_legend_f /
  plot_spacer()+
  plot_layout(heights = c(0.9, 0.025, 0.05, 0.025)) +
  plot_annotation(tag_levels = NULL)

ggsave("./03_plots/copilot/DE_DTU_clock_targets_in_clusters_final_fig.png", grid_2x2_legend, width = 8, height = 12, dpi = 300)

# Rhythmicity----
# *DE----
# Join two data frames
DE_df_d1d2_rhy <- DE_DTU_clock %>%
  filter(type == 'DE') %>% 
  left_join(DE_d1d2_rhy, by = "cluster") %>%
  mutate(dataset = "d1d2")

DE_df_d1d5_rhy <- DE_DTU_clock %>%
  filter(type == 'DE') %>%
  left_join(DE_d1d5_rhy, by = "cluster") %>%
  mutate(dataset = "d1d5")

head(DE_df_d1d5_rhy)

# Count genes per cluster × clock gene
DE_df_d1d2_rhy_counts <- DE_df_d1d2_rhy %>%
  dplyr::count(rhy_category, cluster, clock)

DE_df_d1d5_rhy_counts <- DE_df_d1d5_rhy %>%
  dplyr::count(rhy_category, cluster, clock)

# Compute total counts per cluster (for sorting)
DE_cluster_rhy_totals_d1d2 <- DE_df_d1d2_rhy_counts %>%
  group_by(rhy_category, cluster) %>%
  summarise(total = sum(n), .groups = "drop")

DE_cluster_rhy_totals_d1d5 <- DE_df_d1d5_rhy_counts %>%
  group_by(rhy_category, cluster) %>%
  summarise(total = sum(n), .groups = "drop")

# Merge totals back and sort clusters within each facet
DE_df_d1d2_rhy_plot <- DE_df_d1d2_rhy_counts %>%
  left_join(DE_cluster_rhy_totals_d1d2, by = c("rhy_category", "cluster")) %>%
  mutate(
    clock = factor(clock, levels = clock_order),
    cluster = factor(cluster, levels = DE_cluster_order),
    dataset = 'DE_d1d2'
  )

DE_df_d1d5_rhy_plot <- DE_df_d1d5_rhy_counts %>%
  left_join(DE_cluster_rhy_totals_d1d5, by = c("rhy_category", "cluster")) %>%
  mutate(
    clock = factor(clock, levels = clock_order),
    cluster = factor(cluster, levels = DE_cluster_order),
    dataset = 'DE_d1d5'
  )

DE_df_rhy_all <- bind_rows(DE_df_d1d2_rhy_plot, DE_df_d1d5_rhy_plot) %>% 
  mutate(dataset = case_when(dataset == 'DE_d1d2' ~ 'DE d1d2',
                             dataset == 'DE_d1d5' ~ 'DE d1d5',
                             TRUE ~ dataset))

DE_rhy_stacked_plot <- ggplot(DE_df_rhy_all, aes(x = cluster, y = n, fill = clock)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = tf_palette, breaks = clock_order) +
  scale_x_reordered() +
  facet_grid(rhy_category ~ dataset, scales = "free_x") +
  labs(
    x = "Cluster",
    y = "Gene count",
    fill = "Clock context"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box.spacing = unit(2.5, "pt"),
    legend.box.background = element_rect(color = "black")
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) + 
  labs(tag = "B")

DE_rhy_final_figure <- 
  DE_cluster_size_plot /
  phase_strip /
  DE_rhy_stacked_plot +
  plot_layout(heights = c(0.1, 0.03, 0.87)) 

ggsave("./03_plots/copilot/DE_rhy_clock_targets_in_clusters.png", DE_rhy_final_figure, width = 6, height = 12, dpi = 300)

# *DTU----
# Join two data frames
DTU_df_d1d2_rhy <- DE_DTU_clock %>%
  filter(type == 'DTU') %>% 
  left_join(DTU_d1d2_rhy, by = "cluster") %>%
  mutate(dataset = "d1d2")

DTU_df_d1d5_rhy <- DE_DTU_clock %>%
  filter(type == 'DTU') %>%
  left_join(DTU_d1d5_rhy, by = "cluster") %>%
  mutate(dataset = "d1d5")

head(DTU_df_d1d5_rhy)

# Count genes per cluster × clock gene
DTU_df_d1d2_rhy_counts <- DTU_df_d1d2_rhy %>%
  dplyr::count(rhy_category, cluster, clock)

DTU_df_d1d5_rhy_counts <- DTU_df_d1d5_rhy %>%
  dplyr::count(rhy_category, cluster, clock)

# Compute total counts per cluster (for sorting)
DTU_cluster_rhy_totals_d1d2 <- DTU_df_d1d2_rhy_counts %>%
  group_by(rhy_category, cluster) %>%
  summarise(total = sum(n), .groups = "drop")

DTU_cluster_rhy_totals_d1d5 <- DTU_df_d1d5_rhy_counts %>%
  group_by(rhy_category, cluster) %>%
  summarise(total = sum(n), .groups = "drop")

# Merge totals back and sort clusters within each facet
DTU_df_d1d2_rhy_plot <- DTU_df_d1d2_rhy_counts %>%
  left_join(DTU_cluster_rhy_totals_d1d2, by = c("rhy_category", "cluster")) %>%
  mutate(
    clock = factor(clock, levels = clock_order),
    cluster = factor(cluster, levels = DTU_cluster_order),
    dataset = 'DTU_d1d2'
  )

DTU_df_d1d5_rhy_plot <- DTU_df_d1d5_rhy_counts %>%
  left_join(DTU_cluster_rhy_totals_d1d5, by = c("rhy_category", "cluster")) %>%
  mutate(
    clock = factor(clock, levels = clock_order),
    cluster = factor(cluster, levels = DTU_cluster_order),
    dataset = 'DTU_d1d5'
  )

DTU_df_rhy_all <- bind_rows(DTU_df_d1d2_rhy_plot, DTU_df_d1d5_rhy_plot) %>% 
  mutate(dataset = case_when(dataset == 'DTU_d1d2' ~ 'DTU d1d2',
                             dataset == 'DTU_d1d5' ~ 'DTU d1d5',
                             TRUE ~ dataset))

DTU_rhy_stacked_plot <- ggplot(DTU_df_rhy_all, aes(x = cluster, y = n, fill = clock)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = tf_palette, breaks = clock_order) +
  scale_x_reordered() +
  facet_grid(rhy_category ~ dataset, scales = "free_x") +
  labs(
    x = "Cluster",
    y = "Gene count",
    fill = "Clock context"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.box.spacing = unit(2.5, "pt"),
    legend.box.background = element_rect(color = "black")
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) + 
  labs(tag = "B")

DTU_rhy_final_figure <- 
  DTU_cluster_size_plot /
  phase_strip /
  DTU_rhy_stacked_plot +
  plot_layout(heights = c(0.1, 0.03, 0.87)) 

ggsave("./03_plots/copilot/DTU_rhy_clock_targets_in_clusters.png", DTU_rhy_final_figure, width = 6, height = 12, dpi = 300)

# patched plot - amp and rhy----
DE_stacked_noleg <- DE_stacked_plot + theme(legend.position = "none")
DTU_stacked_noleg <- DTU_stacked_plot + theme(legend.position = "none")
DE_rhy_stacked_noleg <- DE_rhy_stacked_plot + theme(legend.position = "none")
DTU_rhy_stacked_noleg <- DTU_rhy_stacked_plot + theme(legend.position = "none")

shared_legend_f <- extract_legend(DE_stacked_plot)

top_row_amp_rhy <- 
  (DE_cluster_size_plot + labs(tag = "A")) |
  (DTU_cluster_size_plot + labs(tag = "D"))

middle_row_amp_rhy <- 
  ((DE_stacked_noleg + labs(tag = "B")) |
  (DTU_stacked_noleg + labs(tag = "E"))) +
  plot_annotation('Amplitude')

bottom_row_amp_rhy <- 
  ((DE_rhy_stacked_noleg + labs(tag = "C")) |
  (DTU_rhy_stacked_noleg + labs(tag = "F"))) +
  plot_annotation('Rhythmicity')

grid_2x3 <- 
  top_row_amp_rhy /
  phase_strip /
  middle_row_amp_rhy +
  bottom_row_amp_rhy +
  plot_layout(widths = c(1,1), heights = c(0.2, 0.03, 1, 1))

grid_2x3_legend <- 
  wrap_plots(grid_2x3) /
  plot_spacer() /
  shared_legend_f /
  plot_spacer()+
  plot_layout(heights = c(0.9, 0.025, 0.05, 0.025)) +
  plot_annotation(tag_levels = NULL)

ggsave("./03_plots/copilot/DE_DTU_clock_targets_in_clusters_final_fig_amp_rhy.png", grid_2x3_legend, width = 8, height = 16, dpi = 300)


# ------------------------------------------------
# Enrichment – Odds Ratios – Heatmap 
# ------------------------------------------------
# Amplitude----
# *DE heatmap----
clock_all <- DE_DTU_clock %>%
  mutate(cluster = as.numeric(cluster)) %>%
  mutate(clock = case_when(
    clock %in% c("CCA1-Nagel", "CCA1-Kamioka") ~ "CCA1",
    TRUE ~ clock
  ))

clock_order_hm <- c(
  "CCA1", "LHY",
  "PRR5", "PRR7",
  "TOC1",
  "LUX", "ELF3", "ELF4"
)

compute_enrichment <- function(clock_df, universe_df, type_label) {
  
  universe <- universe_df %>%
    distinct(gene_ID, cluster) %>%
    mutate(cluster = as.numeric(cluster))
  
  targets <- clock_df %>%
    distinct(gene_ID, clock) %>%
    mutate(clock = as.character(clock))
  
  cluster_sizes <- universe %>%
    group_by(cluster) %>%
    summarise(n_in_cluster = n(), .groups = "drop")
  
  clock_sizes <- targets %>%
    group_by(clock) %>%
    summarise(n_targets = n(), .groups = "drop")
  
  a_tbl <- universe %>%
    inner_join(targets, by = "gene_ID") %>%
    group_by(cluster, clock) %>%
    summarise(a = n(), .groups = "drop")
  
  all_combos <- expand.grid(
    cluster = unique(universe$cluster),
    clock   = unique(targets$clock),
    stringsAsFactors = FALSE
  ) %>% as_tibble()
  
  enrich <- all_combos %>%
    left_join(a_tbl,         by = c("cluster", "clock")) %>%
    left_join(cluster_sizes, by = "cluster") %>%
    left_join(clock_sizes,   by = "clock") %>%
    mutate(
      a = ifelse(is.na(a), 0, a),
      n_in_cluster = ifelse(is.na(n_in_cluster), 0, n_in_cluster),
      n_targets    = ifelse(is.na(n_targets), 0, n_targets),
      b = n_in_cluster - a,
      c = n_targets - a,
      d = nrow(universe) - a - b - c
    ) %>%
    rowwise() %>%
    mutate(
      p = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
      OR = (a * d) / ((b + 1e-9) * (c + 1e-9)),
      log2_OR = log2(OR + 1e-9),
      type = type_label
    ) %>%
    ungroup()
  
  enrich
}

# Step 3 — enrichment with unified cluster map
enrich_DE <- compute_enrichment(
  clock_df    = clock_all %>% filter(type == "DE"),
  universe_df = calixto_S7A %>% select(gene_ID, cluster),
  type_label  = "DE"
) %>%
  group_by(clock) %>%
  mutate(FDR = p.adjust(p, "BH")) %>%
  ungroup()

# Step 4 — log2(OR) matrix
mat_DE <- enrich_DE %>%
  group_by(cluster, clock) %>%
  summarise(log2_OR = mean(log2_OR), .groups = "drop") %>%
  mutate(
    cluster = as.numeric(cluster),
    clock   = factor(clock, levels = clock_order_hm)
  ) %>%
  tidyr::pivot_wider(names_from = clock, values_from = log2_OR) %>%
  arrange(cluster) %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix()

mat_DE <- mat_DE[, clock_order_hm]
head(mat_DE)

# FDR matrix
FDR_DE <- enrich_DE %>%
  select(cluster, clock, FDR) %>%
  mutate(
    cluster = as.numeric(cluster),
    clock   = factor(clock, levels = clock_order_hm)
  ) %>%
  tidyr::pivot_wider(names_from = clock, values_from = FDR) %>%
  arrange(cluster) %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix()

FDR_DE <- FDR_DE[, clock_order_hm]

# significance mask for outlines (FDR < 0.05 and positive log2_OR)
sig_DE_outline <- (FDR_DE < 0.05) & (mat_DE > 0)

# lock row order
mat_DE         <- mat_DE[order(as.numeric(rownames(mat_DE))), ]
FDR_DE         <- FDR_DE[rownames(mat_DE), ]
sig_DE_outline <- sig_DE_outline[rownames(mat_DE), ]

# amplitude annotations
collapse_amp <- function(df) {
  df %>%
    mutate(cluster = as.numeric(cluster)) %>%
    group_by(cluster) %>%
    summarise(amp = names(sort(table(amp_category), decreasing = TRUE))[1]) %>%
    arrange(cluster)
}

amp_DE_d1d2 <- collapse_amp(DE_df_d1d2)
amp_DE_d1d5 <- collapse_amp(DE_df_d1d5)

amp_vec_DE_d1d2 <- amp_DE_d1d2$amp
amp_vec_DE_d1d5 <- amp_DE_d1d5$amp

amp_cols <- c(
  "gain high"   = "#006400",
  "gain medium" = "#66c2a5",
  "other"       = "#bdbdbd",
  "lose medium" = "#fc8d62",
  "lose high"   = "#8b0000"
)

cluster_sizes_DE <- table(calixto_S7A$cluster)

cluster_sizes_DE_mat <- matrix(
  cluster_sizes_DE[as.character(rownames(mat_DE))],
  ncol = 1,
  dimnames = list(rownames(mat_DE), "size")
)

row_hmap_DE <- rowAnnotation(
  "d1 vs d2" = amp_vec_DE_d1d2,
  "d1 vs d5" = amp_vec_DE_d1d5,
  size       = anno_barplot(
    cluster_sizes_DE_mat,
    gp = gpar(fill = "grey60", col = NA),
    border = FALSE,
    axis_param = list(labels = FALSE)
  ),
  col = list(
    "d1 vs d2" = amp_cols,
    "d1 vs d5" = amp_cols
  ),
  gap = unit(2, "mm"),
  show_legend = FALSE
)


# colour scale
col_fun <- colorRamp2(
  c(0, 1, 2, 3, 4),
  c("white", "#bdd7e7", "#6baed6", "#3182bd", "#08519c")
)

# heatmap
h_DE <- Heatmap(
  mat_DE,
  name = "log2(OR)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = rownames(mat_DE),
  column_order = colnames(mat_DE),
  #column_names_rot = 45,
  right_annotation = row_hmap_DE,
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    at = c(0, 1, 2, 3, 4),
    labels = c("0", "1", "2", "3", "4"),
    color_bar = "continuous"
  ),
  use_raster = FALSE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    
    # outline + number for significant positive tiles
    if (sig_DE_outline[i, j]) {
      grid.rect(x, y, w, h, gp = gpar(fill = NA, col = "black", lwd = 1))
      val     <- mat_DE[i, j]
      txt_col <- if (val > 2.5) "white" else "black"
      grid.text(
        sprintf("%.1f", val),
        x, y,
        gp = gpar(fontsize = 16, col = txt_col)
      )
    }
  }
)

# legends
lgd_list_amp <- list(
  Legend(
    labels = c("Gain High", "Gain Medium", "Other", "Lose Medium", "Lose High"),
    title = "cluster context",
    type = "points",
    pch = 15,
    title_gp = gpar(fontsize = 14),
    size = unit(5, "mm"),
    border = "black",
    legend_gp = gpar(col = amp_cols)
  )
)

draw(h_DE, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_amp)
# dev.off()

# *DTU heatmap----
# 1. Compute DTU enrichment
enrich_DTU <- compute_enrichment(
  clock_df    = clock_all %>% filter(type == "DTU"),
  universe_df = calixto_S7B %>% select(gene_ID, cluster),
  type_label  = "DTU"
) %>%
  group_by(clock) %>%
  mutate(FDR = p.adjust(p, "BH")) %>%
  ungroup()

# 2. Build log₂(OR) matrix for DTU
mat_DTU <- enrich_DTU %>%
  group_by(cluster, clock) %>%
  summarise(log2_OR = mean(log2_OR), .groups = "drop") %>%
  mutate(
    cluster = as.numeric(cluster),
    clock   = factor(clock, levels = clock_order_hm)
  ) %>%
  tidyr::pivot_wider(names_from = clock, values_from = log2_OR) %>%
  arrange(cluster) %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix()

mat_DTU <- mat_DTU[, clock_order_hm]

# 3. Build FDR matrix
FDR_DTU <- enrich_DTU %>%
  select(cluster, clock, FDR) %>%
  mutate(
    cluster = as.numeric(cluster),
    clock   = factor(clock, levels = clock_order_hm)
  ) %>%
  tidyr::pivot_wider(names_from = clock, values_from = FDR) %>%
  arrange(cluster) %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix()

FDR_DTU <- FDR_DTU[, clock_order_hm]

# 4. Significance mask (outline only)
sig_DTU_outline <- (FDR_DTU < 0.05) & (mat_DTU > 0)

# 5. Lock row order
mat_DTU         <- mat_DTU[order(as.numeric(rownames(mat_DTU))), ]
FDR_DTU         <- FDR_DTU[rownames(mat_DTU), ]
sig_DTU_outline <- sig_DTU_outline[rownames(mat_DTU), ]

# 6. Build amplitude annotations (DTU)
amp_DTU_d1d2 <- collapse_amp(DTU_df_d1d2)
amp_DTU_d1d5 <- collapse_amp(DTU_df_d1d5)

amp_vec_DTU_d1d2 <- amp_DTU_d1d2$amp
amp_vec_DTU_d1d5 <- amp_DTU_d1d5$amp

cluster_sizes_DTU <- table(calixto_S7B$cluster)

cluster_sizes_DTU_mat <- matrix(
  cluster_sizes_DTU[as.character(rownames(mat_DTU))],
  ncol = 1,
  dimnames = list(rownames(mat_DTU), "size")
)

row_hmap_DTU <- rowAnnotation(
  "d1 vs d2" = amp_vec_DTU_d1d2,
  "d1 vs d5" = amp_vec_DTU_d1d5,
  size       = anno_barplot(
    cluster_sizes_DTU_mat,
    gp = gpar(fill = "grey60", col = NA),
    border = FALSE,
    axis_param = list(labels = FALSE)
  ),
  col = list(
    "d1 vs d2" = amp_cols,
    "d1 vs d5" = amp_cols
  ),
  gap = unit(2, "mm"),
  show_legend = FALSE
)

# 7. Colour scale (identical to DE)
col_fun <- colorRamp2(
  c(0, 1, 2, 3, 4),
  c("white", "#bdd7e7", "#6baed6", "#3182bd", "#08519c")
)

# 8. DTU Heatmap (identical design)
h_DTU <- Heatmap(
  mat_DTU,
  name = "log2(OR)",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = rownames(mat_DTU),
  column_order = colnames(mat_DTU),
  #column_names_rot = 45,
  right_annotation = row_hmap_DTU,
  
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    at = c(0, 1, 2, 3, 4),
    labels = c("0", "1", "2", "3", "4"),
    color_bar = "continuous"
  ),
  
  use_raster = FALSE,
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    
    if (sig_DTU_outline[i, j]) {
      grid.rect(x, y, w, h, gp = gpar(fill = NA, col = "black", lwd = 1))
      
      val     <- mat_DTU[i, j]
      txt_col <- if (val > 2.5) "white" else "black"
      
      grid.text(
        sprintf("%.1f", val),
        x, y,
        gp = gpar(fontsize = 16, col = txt_col)
      )
    }
  }
)

# # 9. Legend (same as DE)
# lgd_list <- list(
#   Legend(
#     labels = c("Gain High", "Gain Medium", "Other", "Lose Medium", "Lose High"),
#     title = "cluster context",
#     type = "points",
#     pch = 15,
#     title_gp = gpar(fontsize = 10),
#     size = unit(5, "mm"),
#     border = "black",
#     legend_gp = gpar(col = amp_cols)
#   )
# )

draw(h_DTU, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_amp)

# Rhythmicity----
# * DE heatmap----
# rhythmicity annotations
collapse_rhy <- function(df) {
  df %>%
    mutate(cluster = as.numeric(cluster)) %>%
    group_by(cluster) %>%
    summarise(rhy = names(sort(table(rhy_category), decreasing = TRUE))[1]) %>%
    arrange(cluster)
}

rhy_DE_d1d2 <- collapse_rhy(DE_df_d1d2_rhy)
rhy_DE_d1d5 <- collapse_rhy(DE_df_d1d5_rhy)

rhy_vec_DE_d1d2 <- rhy_DE_d1d2$rhy
rhy_vec_DE_d1d5 <- rhy_DE_d1d5$rhy

rhy_cols <- c(
  "nr-r"   = "#66c2a5",
  "r-nr"   = "#fc8d62",
  "r-r"    = "#006400",
  "nr-nr"  = "#bdbdbd"
)

row_hmap_DE_rhy <- rowAnnotation(
  "d1 vs d2" = rhy_vec_DE_d1d2,
  "d1 vs d5" = rhy_vec_DE_d1d5,
  size       = anno_barplot(
    cluster_sizes_DE_mat,
    gp = gpar(fill = "grey60", col = NA),
    border = FALSE,
    axis_param = list(labels = FALSE)
  ),
  col = list(
    "d1 vs d2" = rhy_cols,
    "d1 vs d5" = rhy_cols
  ),
  gap = unit(2, "mm"),
  show_legend = FALSE
)

# heatmap
h_DE_rhy <- Heatmap(
  mat_DE,
  name = "log2(OR)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = rownames(mat_DE),
  column_order = colnames(mat_DE),
  #column_names_rot = 45,
  right_annotation = row_hmap_DE_rhy,
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    at = c(0, 1, 2, 3, 4),
    labels = c("0", "1", "2", "3", "4"),
    color_bar = "continuous"
  ),
  use_raster = FALSE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    
    # outline + number for significant positive tiles
    if (sig_DE_outline[i, j]) {
      grid.rect(x, y, w, h, gp = gpar(fill = NA, col = "black", lwd = 1))
      val     <- mat_DE[i, j]
      txt_col <- if (val > 2.5) "white" else "black"
      grid.text(
        sprintf("%.1f", val),
        x, y,
        gp = gpar(fontsize = 16, col = txt_col)
      )
    }
  }
)

# legends
lgd_list_rhy <- list(
  Legend(
    labels = c("nr-r", "r-nr", "r-r", "nr-nr"),
    title = "cluster context",
    type = "points",
    pch = 15,
    title_gp = gpar(fontsize = 14),
    size = unit(5, "mm"),
    border = "black",
    legend_gp = gpar(col = rhy_cols)
  )
)

draw(h_DE_rhy, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_rhy)

# * DTU heatmap----
# rhythmicity annotations

rhy_DTU_d1d2 <- collapse_rhy(DTU_df_d1d2_rhy)
rhy_DTU_d1d5 <- collapse_rhy(DTU_df_d1d5_rhy)

rhy_vec_DTU_d1d2 <- rhy_DTU_d1d2$rhy
rhy_vec_DTU_d1d5 <- rhy_DTU_d1d5$rhy

row_hmap_DTU_rhy <- rowAnnotation(
  "d1 vs d2" = rhy_vec_DTU_d1d2,
  "d1 vs d5" = rhy_vec_DTU_d1d5,
  size       = anno_barplot(
    cluster_sizes_DTU_mat,
    gp = gpar(fill = "grey60", col = NA),
    border = FALSE,
    axis_param = list(labels = FALSE)
  ),
  col = list(
    "d1 vs d2" = rhy_cols,
    "d1 vs d5" = rhy_cols
  ),
  gap = unit(2, "mm"),
  show_legend = FALSE
)

# heatmap
h_DTU_rhy <- Heatmap(
  mat_DTU,
  name = "log2(OR)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = rownames(mat_DTU),
  column_order = colnames(mat_DTU),
  #column_names_rot = 45,
  right_annotation = row_hmap_DTU_rhy,
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    at = c(0, 1, 2, 3, 4),
    labels = c("0", "1", "2", "3", "4"),
    color_bar = "continuous"
  ),
  use_raster = FALSE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    
    # outline + number for significant positive tiles
    if (sig_DTU_outline[i, j]) {
      grid.rect(x, y, w, h, gp = gpar(fill = NA, col = "black", lwd = 1))
      val     <- mat_DTU[i, j]
      txt_col <- if (val > 2.5) "white" else "black"
      grid.text(
        sprintf("%.1f", val),
        x, y,
        gp = gpar(fontsize = 16, col = txt_col)
      )
    }
  }
)

draw(h_DTU_rhy, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_rhy)

# png version
# png("./03_plots/h_DE.png", width = 2000, height = 2500, res = 300)
# draw(h_DE, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_amp)
# dev.off()
# 
# png("./03_plots/h_DTU.png", width = 2000, height = 2500, res = 300)
# draw(h_DTU, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_amp)
# dev.off()
# 
# png("./03_plots/h_DE_rhy.png", width = 2000, height = 2500, res = 300)
# draw(h_DE_rhy, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_rhy)
# dev.off()
# 
# png("./03_plots/h_DTU_rhy.png", width = 2000, height = 2500, res = 300)
# draw(h_DTU_rhy, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_rhy)
# dev.off()

# pdf version
pdf("./03_plots/h_DE.pdf", width = 7, height = 9)
draw(h_DE, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_amp)
dev.off()

pdf("./03_plots/h_DTU.pdf", width = 7, height = 9)
draw(h_DTU, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_amp)
dev.off()

pdf("./03_plots/h_DE_rhy.pdf", width = 7, height = 9)
draw(h_DE_rhy, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_rhy)
dev.off()

pdf("./03_plots/h_DTU_rhy.pdf", width = 7, height = 9)
draw(h_DTU_rhy, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list_rhy)
dev.off()

# read back in png version
# img_DE      <- image_read("./03_plots/h_DE.png")
# img_DTU     <- image_read("./03_plots/h_DTU.png")
# img_DE_rhy  <- image_read("./03_plots/h_DE_rhy.png")
# img_DTU_rhy <- image_read("./03_plots/h_DTU_rhy.png")

# read back in pdf version
img_DE      <- magick::image_read_pdf("./03_plots/h_DE.pdf")
img_DTU     <- magick::image_read_pdf("./03_plots/h_DTU.pdf")
img_DE_rhy  <- magick::image_read_pdf("./03_plots/h_DE_rhy.pdf")
img_DTU_rhy <- magick::image_read_pdf("./03_plots/h_DTU_rhy.pdf")

pad_image <- function(img, pad = 120) {
  w <- magick::image_info(img)$width
  h <- magick::image_info(img)$height
  
  # Create a padded canvas
  canvas <- magick::image_blank(
    width  = w + 2*pad,
    height = h + 2*pad,
    color  = "white"
  )
  
  # Composite heatmap onto the padded canvas
  magick::image_composite(canvas, img, offset = paste0("+", pad, "+", pad))
}

img_DE      <- pad_image(img_DE)
img_DTU     <- pad_image(img_DTU)
img_DE_rhy  <- pad_image(img_DE_rhy)
img_DTU_rhy <- pad_image(img_DTU_rhy)


label_panel <- function(img, label) {
  magick::image_annotate(
    img,
    label,
    size = 90,
    weight = 700,
    gravity = "northwest",
    location = "+20+20"   # padding from top-left
  )
}

img_DE      <- label_panel(img_DE, "A")
img_DTU     <- label_panel(img_DTU, "B")
img_DE_rhy  <- label_panel(img_DE_rhy, "C")
img_DTU_rhy <- label_panel(img_DTU_rhy, "D")

make_row_label <- function(text, width) {
  magick::image_blank(width = width, height = 120, color = "white") %>%
    magick::image_annotate(
      text,
      size = 90,
      weight = 700,
      gravity = "west",
      location = "+40+0"   # left padding
    )
}

row_width <- max(magick::image_info(img_DE)$width + magick::image_info(img_DTU)$width)

amp_label <- make_row_label("Amplitude", row_width)
rhy_label <- make_row_label("Rhythmicity", row_width)

top_row    <- magick::image_append(c(img_DE, img_DTU))
bottom_row <- magick::image_append(c(img_DE_rhy, img_DTU_rhy))

final_4panel <- magick::image_append(
  c(amp_label, top_row, rhy_label, bottom_row),
  stack = TRUE
)

# png version
#image_write(final_4panel, "./03_plots/final_4panel_heatmaps.png")

# pdf version
#magick::image_write(final_4panel, "./03_plots/final_4panel_heatmaps.pdf")

magick::image_write(final_4panel, path = "./03_plots/final_4panel_heatmaps.pdf", format = "pdf")

# 2-panel Amplitude heatmap _ DE and DTU----

# Enlarged cluster‑context legend
lgd_list_big <- list(
  Legend(
    labels = c("Gain High", "Gain Medium", "Other", "Lose Medium", "Lose High"),
    title = "cluster context",
    type = "points",
    pch = 15,
    ncol=5,
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    size = unit(5, "mm"),
    #border = "black",
    legend_gp = gpar(col = amp_cols),
    row_gap = unit(8, "mm")
  )
)

# Draw each amplitude heatmap to PDF (vector, crisp)
#pdf("./03_plots/amp_DE.pdf", width = 7, height = 9)
pdf("./03_plots/amp_DE.pdf", width = 5, height = 5)
draw(
  h_DE,
  heatmap_legend_side = "bottom",
  annotation_legend_list = lgd_list_big,
  annotation_legend_side = "bottom")
# draw(h_DE, 
#      heatmap_legend_side = "bottom")
dev.off()

#pdf("./03_plots/amp_DTU.pdf", width = 7, height = 9)
pdf("./03_plots/amp_DTU.pdf", width = 5, height = 5)
draw(
  h_DTU,
  heatmap_legend_side = "bottom",
  annotation_legend_list = lgd_list_big,
  annotation_legend_side = "bottom")
# draw(h_DTU,
#      heatmap_legend_side = "bottom")
dev.off()

# Read PDFs back in
img_DE  <- magick::image_read_pdf("./03_plots/amp_DE.pdf")
img_DTU <- magick::image_read_pdf("./03_plots/amp_DTU.pdf")

# Add padding around each heatmap
pad_image <- function(img, pad = 120) {
  w <- magick::image_info(img)$width
  h <- magick::image_info(img)$height
  
  canvas <- magick::image_blank(
    width  = w + 2*pad,
    height = h + 2*pad,
    color  = "white"
  )
  
  magick::image_composite(canvas, img, offset = paste0("+", pad, "+", pad))
}

img_DE  <- pad_image(img_DE)
img_DTU <- pad_image(img_DTU)

# Equalise tile sizes: DTU has 10 clusters, DE has 12
# scale_factor <- 0.85
# new_height <- round(magick::image_info(img_DTU)$height * scale_factor)
# img_DTU <- magick::image_resize(img_DTU, geometry = paste0("x", new_height))


# Add panel labels A and B (bold, 90 pt, left‑aligned)
label_panel <- function(img, label) {
  magick::image_annotate(
    img,
    label,
    size = 90,
    weight = 700,
    gravity = "northwest",
    location = "+20+20"
  )
}

img_DE  <- label_panel(img_DE, "A")
img_DTU <- label_panel(img_DTU, "B")

# Assemble the 2‑panel amplitude figure
amp_2panel <- magick::image_append(c(img_DE, img_DTU), stack = TRUE)

# Save as vector‑quality PDF (for Word)
magick::image_write(amp_2panel,
                    "./03_plots/heatmap_amplitude_main_figure_alt.pdf",
                    format = "pdf")

# Generate the legend as a standalone PDF using a tiny dummy heatmap, 
# then treat it as a separate image panel
# dummy <- Heatmap(
#   matrix(0, 1, 1),
#   name = "dummy",
#   col = c("white" = "white"),
#   show_heatmap_legend = FALSE,
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   show_row_names = FALSE,
#   show_column_names = FALSE,
#   rect_gp = gpar(col = NA, fill = NA)   # ← THIS is the missing piece
# )
# 
# pdf("./03_plots/amp_legend.pdf", width = 5, height = 1.2)
# draw(
#   dummy,
#   annotation_legend_list = lgd_list_big,
#   annotation_legend_side = "bottom",
#   merge_legend = TRUE
# )
# dev.off()

# pdf("./03_plots/amp_legend.pdf", width = 4, height = 1.2)
# draw(dummy,
#      annotation_legend_list = lgd_list_big,
#      heatmap_legend_side = "bottom",
#      merge_legend = TRUE)
# dev.off()
# 
# legend_img <- magick::image_read_pdf("./03_plots/amp_legend.pdf")
# 
# final_amp_2panel <- magick::image_append(
#   c(amp_2panel, legend_img),
#   stack = TRUE
# )

# Save as vector‑quality PDF (for Word)
# magick::image_write(final_amp_2panel,
#                     "./03_plots/heatmap_amplitude_main_figure_legend.pdf",
#                     format = "pdf")


# UNIFIED PERMUTATION + PHASE MODEL WORKFLOW----
# For DE and DTU amplitude gain/loss contrasts
# Using linear and circular phase models
# With per‑TF deltas as the biological unit of replication

# *Setup: TF phases and cluster definitions----
###############################################
# 0. Setup: TF phases and cluster definitions #
###############################################

# ZT phases (hours) from literature
TF_phase_ZT <- c(
  CCA1 = 1,
  LHY  = 1,
  PRR5 = 10,
  PRR7 = 6,
  TOC1 = 12,
  LUX  = 14,
  ELF3 = 14,
  ELF4 = 14
)

# Convert to radians for circular model
TF_theta <- 2 * pi * TF_phase_ZT / 24

cluster_context_DE <- amp_vec_DE_d1d2
names(cluster_context_DE) <- as.character(1:12)

cluster_context_DTU <- amp_vec_DTU_d1d2
names(cluster_context_DTU) <- as.character(1:10)

# Normalise context labels
cluster_context_DE  <- trimws(tolower(cluster_context_DE))
cluster_context_DTU <- trimws(tolower(cluster_context_DTU))

# Gain/loss cluster definitions
gain_clusters <- c("gain high", "gain medium")
loss_clusters <- c("lose medium", "lose high")
loss_clusters_DTU <- "lose medium"

# *Convert cluster‑level log2OR matrices → TF × context matrices----
###############################################
# 1. Convert cluster matrix → TF × context    #
###############################################

make_TF_context_matrix <- function(mat, context_vec) {
  df <- as.data.frame(mat)
  rownames(df) <- names(context_vec)
  
  df$context <- trimws(tolower(context_vec[rownames(df)]))
  
  df %>%
    mutate(cluster = rownames(.)) %>%
    pivot_longer(cols = all_of(colnames(mat)),
                 names_to = "TF", values_to = "log2OR") %>%
    select(TF, context, log2OR) %>%
    pivot_wider(names_from = context, values_from = log2OR,
                values_fn = mean) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("TF")
}

log2OR_mat_DE  <- make_TF_context_matrix(mat_DE,  cluster_context_DE)
log2OR_mat_DTU <- make_TF_context_matrix(mat_DTU, cluster_context_DTU)

colnames(log2OR_mat_DE)
colnames(log2OR_mat_DTU)

# *Compute per‑TF deltas (gain – loss)----
###############################################
# 2. Per‑TF delta calculation                 #
###############################################

compute_TF_deltas <- function(mat, gain_clusters, loss_clusters) {
  gain_vals <- rowMeans(as.matrix(mat[, gain_clusters]), na.rm = TRUE)
  loss_vals <- rowMeans(as.matrix(mat[, loss_clusters]), na.rm = TRUE)
  gain_vals - loss_vals
}

TF_deltas_DE_d1d2  <- compute_TF_deltas(log2OR_mat_DE,  gain_clusters, loss_clusters)
TF_deltas_DTU_d1d2 <- compute_TF_deltas(log2OR_mat_DTU, gain_clusters, loss_clusters_DTU)

# *Linear regression on phase (monotonic trend)----
###############################################
# 3. Linear regression on phase               #
###############################################

run_linear_model <- function(TF_deltas) {
  df <- data.frame(
    delta = as.numeric(TF_deltas),
    phase = TF_phase_ZT[names(TF_deltas)]
  )
  lm(delta ~ phase, data = df)
}

fit_DE_linear_d1d2  <- run_linear_model(TF_deltas_DE_d1d2)
fit_DTU_linear_d1d2 <- run_linear_model(TF_deltas_DTU_d1d2)

# *Permutation test for linear slope----
###############################################
# 4. Permutation test for linear slope        #
###############################################

permute_linear <- function(TF_deltas, B = 10000) {
  df <- data.frame(
    delta = as.numeric(TF_deltas),
    phase = TF_phase_ZT[names(TF_deltas)]
  )
  beta_obs <- coef(lm(delta ~ phase, data = df))["phase"]
  
  beta_perm <- numeric(B)
  for (b in seq_len(B)) {
    df$phase_perm <- sample(df$phase)
    beta_perm[b] <- coef(lm(delta ~ phase_perm, data = df))["phase_perm"]
  }
  
  p_one <- (sum(beta_perm >= beta_obs) + 1) / (B + 1)
  p_two <- (sum(abs(beta_perm) >= abs(beta_obs)) + 1) / (B + 1)
  
  list(beta_obs = beta_obs, p_one = p_one, p_two = p_two)
}

perm_DE_linear  <- permute_linear(TF_deltas_DE_d1d2)
perm_DTU_linear <- permute_linear(TF_deltas_DTU_d1d2)

# *Circular regression (sinusoidal phase model)----
###############################################
# 5. Circular regression (sinusoidal model)   #
###############################################

run_circular_model <- function(TF_deltas) {
  df <- data.frame(
    delta = as.numeric(TF_deltas),
    theta = TF_theta[names(TF_deltas)]
  )
  df$cos_theta <- cos(df$theta)
  df$sin_theta <- sin(df$theta)
  
  fit <- lm(delta ~ cos_theta + sin_theta, data = df)
  
  b1 <- coef(fit)["cos_theta"]
  b2 <- coef(fit)["sin_theta"]
  
  A   <- sqrt(b1^2 + b2^2)                 # amplitude
  phi <- atan2(b2, b1)                     # preferred phase (radians)
  ZT  <- (phi * 24 / (2 * pi)) %% 24       # convert to ZT hours
  
  list(fit = fit, amplitude = A, ZT = ZT)
}

circ_DE  <- run_circular_model(TF_deltas_DE_d1d2)
circ_DTU <- run_circular_model(TF_deltas_DTU_d1d2)

# *Permutation test for circular amplitude----
###############################################
# 6. Permutation test for circular amplitude  #
###############################################

permute_circular <- function(TF_deltas, B = 10000) {
  df <- data.frame(
    delta = as.numeric(TF_deltas),
    theta = TF_theta[names(TF_deltas)]
  )
  
  df$cos_theta <- cos(df$theta)
  df$sin_theta <- sin(df$theta)
  
  fit <- lm(delta ~ cos_theta + sin_theta, data = df)
  b1  <- coef(fit)["cos_theta"]
  b2  <- coef(fit)["sin_theta"]
  A_obs <- sqrt(b1^2 + b2^2)
  
  A_perm <- numeric(B)
  for (b in seq_len(B)) {
    df$theta_perm <- sample(df$theta)
    df$cos_perm   <- cos(df$theta_perm)
    df$sin_perm   <- sin(df$theta_perm)
    
    fit_perm <- lm(delta ~ cos_perm + sin_perm, data = df)
    b1p <- coef(fit_perm)["cos_perm"]
    b2p <- coef(fit_perm)["sin_perm"]
    A_perm[b] <- sqrt(b1p^2 + b2p^2)
  }
  
  p <- (sum(A_perm >= A_obs) + 1) / (B + 1)
  list(A_obs = A_obs, p = p)
}

perm_DE_circ  <- permute_circular(TF_deltas_DE_d1d2)
perm_DTU_circ <- permute_circular(TF_deltas_DTU_d1d2)

# streamlined process as functions for DE and DTU x d1d2 and d1d5----
cluster_context_DE_d1d2 <- amp_vec_DE_d1d2
names(cluster_context_DE_d1d2) <- as.character(1:12)

cluster_context_DE_d1d5 <- amp_vec_DE_d1d5
names(cluster_context_DE_d1d5) <- as.character(1:12)

cluster_context_DTU_d1d2 <- amp_vec_DTU_d1d2
names(cluster_context_DTU_d1d2) <- as.character(1:10)

cluster_context_DTU_d1d5 <- amp_vec_DTU_d1d5
names(cluster_context_DTU_d1d5) <- as.character(1:10)

cluster_context_DE_d1d2  <- trimws(tolower(cluster_context_DE_d1d2))
cluster_context_DE_d1d5  <- trimws(tolower(cluster_context_DE_d1d5))

cluster_context_DTU_d1d2 <- trimws(tolower(cluster_context_DTU_d1d2))
cluster_context_DTU_d1d5 <- trimws(tolower(cluster_context_DTU_d1d5))

log2OR_mat_DE_d1d2  <- make_TF_context_matrix(mat_DE,  cluster_context_DE_d1d2)
log2OR_mat_DTU_d1d2 <- make_TF_context_matrix(mat_DTU, cluster_context_DTU_d1d2)

log2OR_mat_DE_d1d5  <- make_TF_context_matrix(mat_DE,  cluster_context_DE_d1d5)
log2OR_mat_DTU_d1d5 <- make_TF_context_matrix(mat_DTU, cluster_context_DTU_d1d5)

compute_TF_deltas_flexible <- function(mat, gain_candidates, loss_candidates) {
  
  available <- colnames(mat)
  
  gain_use <- intersect(gain_candidates, available)
  loss_use <- intersect(loss_candidates, available)
  
  if (length(loss_use) == 0)
    stop("No loss clusters found in matrix.")
  
  if (length(gain_use) == 0) {
    message("No gain clusters found — using gain = 0 for all TFs.")
    gain_vals <- rep(0, nrow(mat))
  } else {
    gain_vals <- rowMeans(as.matrix(mat[, gain_use, drop=FALSE]), na.rm=TRUE)
  }
  
  loss_vals <- rowMeans(as.matrix(mat[, loss_use, drop=FALSE]), na.rm=TRUE)
  
  gain_vals - loss_vals
}

run_phase_analysis <- function(mat, 
                               gain_candidates = c("gain high", "gain medium"),
                               loss_candidates = c("lose high", "lose medium")) {
  
  # 1. Compute deltas
  deltas <- compute_TF_deltas_flexible(mat, gain_candidates, loss_candidates)
  
  # 2. Linear model
  df_lin <- data.frame(
    delta = as.numeric(deltas),
    phase = TF_phase_ZT[names(deltas)]
  )
  fit_lin <- lm(delta ~ phase, data=df_lin)
  beta_obs <- coef(fit_lin)["phase"]
  
  # Permutation
  B <- 10000
  beta_perm <- numeric(B)
  for (b in seq_len(B)) {
    df_lin$phase_perm <- sample(df_lin$phase)
    beta_perm[b] <- coef(lm(delta ~ phase_perm, data=df_lin))["phase_perm"]
  }
  p_lin_one <- (sum(beta_perm >= beta_obs) + 1) / (B + 1)
  p_lin_two <- (sum(abs(beta_perm) >= abs(beta_obs)) + 1) / (B + 1)
  
  # 3. Circular model
  df_circ <- data.frame(
    delta = as.numeric(deltas),
    theta = TF_theta[names(deltas)]
  )
  df_circ$cos_theta <- cos(df_circ$theta)
  df_circ$sin_theta <- sin(df_circ$theta)
  
  fit_circ <- lm(delta ~ cos_theta + sin_theta, data=df_circ)
  b1 <- coef(fit_circ)["cos_theta"]
  b2 <- coef(fit_circ)["sin_theta"]
  A_obs <- sqrt(b1^2 + b2^2)
  ZT_pref <- (atan2(b2, b1) * 24 / (2*pi)) %% 24
  
  # Circular permutation
  A_perm <- numeric(B)
  for (b in seq_len(B)) {
    df_circ$theta_perm <- sample(df_circ$theta)
    df_circ$cos_perm <- cos(df_circ$theta_perm)
    df_circ$sin_perm <- sin(df_circ$theta_perm)
    fit_perm <- lm(delta ~ cos_perm + sin_perm, data=df_circ)
    b1p <- coef(fit_perm)["cos_perm"]
    b2p <- coef(fit_perm)["sin_perm"]
    A_perm[b] <- sqrt(b1p^2 + b2p^2)
  }
  p_circ <- (sum(A_perm >= A_obs) + 1) / (B + 1)
  
  list(
    deltas = deltas,
    linear = list(
      fit = fit_lin,
      beta = beta_obs,
      p_one = p_lin_one,
      p_two = p_lin_two
    ),
    circular = list(
      fit = fit_circ,
      amplitude = A_obs,
      preferred_ZT = ZT_pref,
      p = p_circ
    )
  )
}

res_DE_d1d2 <- run_phase_analysis(
  log2OR_mat_DE_d1d2,
  gain_candidates = c("gain high", "gain medium"),
  loss_candidates = c("lose high", "lose medium")
)

res_DE_d1d5 <- run_phase_analysis(
  log2OR_mat_DE_d1d5,
  gain_candidates = c("gain high", "gain medium"),  # none will be found
  loss_candidates = c("lose high", "lose medium")
)

res_DTU_d1d2 <- run_phase_analysis(
  log2OR_mat_DTU_d1d2,
  gain_candidates = c("gain high", "gain medium"),
  loss_candidates = c("lose high", "lose medium")
)

res_DTU_d1d5 <- run_phase_analysis(
  log2OR_mat_DTU_d1d5,
  gain_candidates = c("gain high", "gain medium"),
  loss_candidates = c("lose high", "lose medium")
)

# * summarise results----
# Build a DE vs DTU × d1d2 vs d1d5 summary table
summarise_phase_result <- function(res, label_condition, label_type) {
  data.frame(
    condition   = label_condition,          # "d1d2" or "d1d5"
    type        = label_type,               # "DE" or "DTU"
    beta_linear = unname(res$linear$beta),
    p_lin_one   = res$linear$p_one,
    p_lin_two   = res$linear$p_two,
    A_circ      = res$circular$amplitude,
    ZT_pref     = res$circular$preferred_ZT,
    p_circ      = res$circular$p
  )
}

summary_all <- rbind(
  summarise_phase_result(res_DE_d1d2,  "d1d2", "DE"),
  summarise_phase_result(res_DTU_d1d2, "d1d2", "DTU"),
  summarise_phase_result(res_DE_d1d5,  "d1d5", "DE"),
  summarise_phase_result(res_DTU_d1d5, "d1d5", "DTU")
)

summary_all

# *Phase–trend plots for all conditions----
make_phase_df <- function(res, label_condition, label_type) {
  data.frame(
    TF      = names(res$deltas),
    delta   = as.numeric(res$deltas),
    phase   = TF_phase_ZT[names(res$deltas)],
    condition = label_condition,
    type      = label_type
  )
}

df_plot <- rbind(
  make_phase_df(res_DE_d1d2,  "d1d2", "DE"),
  make_phase_df(res_DTU_d1d2, "d1d2", "DTU"),
  make_phase_df(res_DE_d1d5,  "d1d5", "DE"),
  make_phase_df(res_DTU_d1d5, "d1d5", "DTU")
)

# Base R version (4 panels)
par(mfrow = c(2, 2))

for (cond in c("d1d2", "d1d5")) {
  for (typ in c("DE", "DTU")) {
    sub <- subset(df_plot, condition == cond & type == typ)
    fit <- lm(delta ~ phase, data = sub)
    
    plot(sub$phase, sub$delta,
         main = paste(typ, cond),
         xlab = "Phase (ZT hours)",
         ylab = "Gain - Loss (delta)",
         pch = 16)
    abline(fit, col = "red", lwd = 2)
    text(sub$phase, sub$delta, labels = sub$TF, pos = 3, cex = 0.8)
  }
}
par(mfrow = c(1, 1))

# Build a unified plotting dataframe----
make_phase_df <- function(res, condition, type) {
  data.frame(
    TF        = names(res$deltas),
    delta     = as.numeric(res$deltas),
    phase     = TF_phase_ZT[names(res$deltas)],
    condition = condition,
    type      = type,
    beta      = res$linear$beta,
    p_lin     = res$linear$p_one,
    A_circ    = res$circular$amplitude,
    ZT_pref   = res$circular$preferred_ZT,
    p_circ    = res$circular$p
  )
}

df_plot <- bind_rows(
  make_phase_df(res_DE_d1d2,  "d1d2", "DE"),
  make_phase_df(res_DTU_d1d2, "d1d2", "DTU"),
  make_phase_df(res_DE_d1d5,  "d1d5", "DE"),
  make_phase_df(res_DTU_d1d5, "d1d5", "DTU")
)

# Create a stats label for each panel
df_stats <- df_plot %>%
  group_by(condition, type) %>%
  summarise(
    beta    = unique(beta),
    p_lin   = unique(p_lin),
    A_circ  = unique(A_circ),
    ZT_pref = unique(ZT_pref),
    p_circ  = unique(p_circ)
  ) %>%
  mutate(
    label = sprintf("β = %.2f\np = %.3f\nZTpref = %.1f",
                    beta, p_lin, ZT_pref)
  )

# The 2×2 ggplot figure
p <- ggplot(df_plot, aes(x = phase, y = delta)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 1) +
  geom_text_repel(aes(label = TF), size = 3.5, max.overlaps = Inf) +
  facet_grid(type ~ condition) +
  theme_bw(base_size = 14) +
  labs(
    x = "TF Phase (ZT hours)",
    y = "Gain – Loss (delta)",
    title = "Phase-dependent amplitude bias across DE/DTU and d1d2/d1d5"
  )

# Add stats to each facet
p_final <- p + geom_text(
  data = df_stats,
  aes(x = -Inf, y = Inf, label = label),
  hjust = -0.1, vjust = 1.1,
  size = 4,
  colour = "black"
)

p_final

ggsave("./03_plots/copilot/phase_trends_2x2.png", p_final, width = 10, height = 8, dpi = 600)


#UpSet plot DE-amp----
kamioka_subset <- DE_clock %>% filter(clock == 'CCA1-Kamioka') %>% mutate(clock = case_when(clock == 'CCA1-Kamioka' ~ 'CCA1', TRUE ~ clock))
nagel_subset <- DE_clock %>% filter(clock == 'CCA1-Nagel') %>% mutate(clock = case_when(clock == 'CCA1-Nagel' ~ 'CCA1', TRUE ~ clock))

cca1_intersection <- inner_join(
  nagel_subset %>% select(gene_ID),
  kamioka_subset %>% select(gene_ID),
  by = "gene_ID"
) %>%
  mutate(clock = "CCA1")

# Convert each ChIP dataset to a vector 
LHY_upset  <- DE_clock %>% filter(clock == 'LHY')  %>% pull(gene_ID)
CCA1_upset <- cca1_intersection %>% pull(gene_ID)
TOC1_upset <- DE_clock %>% filter(clock == 'TOC1') %>% pull(gene_ID)
PRR5_upset <- DE_clock %>% filter(clock == 'PRR5') %>% pull(gene_ID)
PRR7_upset <- DE_clock %>% filter(clock == 'PRR7') %>% pull(gene_ID)
LUX_upset  <- DE_clock %>% filter(clock == 'LUX')  %>% pull(gene_ID)
ELF3_upset <- DE_clock %>% filter(clock == 'ELF3') %>% pull(gene_ID)
ELF4_upset <- DE_clock %>% filter(clock == 'ELF4') %>% pull(gene_ID)

# Build a long-format data frame for ComplexUpset
myGeneSets <- list(LHY = LHY_upset,
                   CCA1 = CCA1_upset,
                   TOC1 = TOC1_upset,
                   PRR5 = PRR5_upset,
                   PRR7 = PRR7_upset,
                   LUX = LUX_upset,
                   ELF3 = ELF3_upset,
                   ELF4 = ELF4_upset)

# Convert to long format
long_df <- enframe(myGeneSets, name = "clock", value = "gene_ID") %>%
  unnest(gene_ID)

# Convert long → wide binary membership matrix
binary_df <- long_df %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = clock,
    values_from = value,
    values_fill = 0
  )

# Build an amplitude metadata table
amp_meta <- DE_d1d2_clock %>%
  select(gene_ID, amp_category) %>%
  distinct() %>%
  rename(amp_d1d2 = amp_category) %>%
  full_join(
    DE_d1d5_clock %>%
      select(gene_ID, amp_category) %>%
      distinct() %>%
      rename(amp_d1d5 = amp_category),
    by = "gene_ID"
  )

# Join amplitude metadata onto the binary_df
binary_df2 <- binary_df %>%
  left_join(amp_meta, by = "gene_ID")

# Fix the factor order (this controls stack order)
binary_df2 <- binary_df2 %>%
  mutate(
    amp_d1d2 = factor(amp_d1d2,
                      levels = c("gain high", "gain medium", "other", "lose medium", "lose high")
    ),
    amp_d1d5 = factor(amp_d1d5,
                      levels = c("gain high", "gain medium", "other", "lose medium", "lose high")
    )
  )

# Define the set names once
# Define the set order exactly as ComplexUpset sees it
set_order <- c("LHY","CCA1","TOC1","PRR5","PRR7","LUX","ELF3","ELF4")

# Now generate numeric intersection codes
binary_df2 <- binary_df2 |>
  dplyr::rowwise() |>
  dplyr::mutate(
    .intersection = {
      present <- c_across(all_of(set_order)) == 1
      paste(which(present), collapse = "-")
    }
  ) |>
  dplyr::ungroup()

# Build the UpSet plot FIRST (this defines the universe)
up <- ComplexUpset::upset(
  binary_df2,
  intersect = c("LHY","CCA1","TOC1","PRR5","PRR7","LUX","ELF3","ELF4"),
  min_size = 5,
  base_annotations = list(
    "Intersection size" = intersection_size()
  ),
  set_sizes=(
    upset_set_size()
    + theme(axis.text.x=element_text(angle=90))
  ),
  queries = list(
    upset_query(set = 'LHY',  fill = '#5aae61'),
    upset_query(set = 'CCA1', fill = '#006d2c'),
    upset_query(set = 'TOC1', fill = '#4d4d4d'),
    upset_query(set = 'PRR5', fill = '#fdae61'),
    upset_query(set = 'PRR7', fill = '#f46d43'),
    upset_query(set = 'LUX',  fill = '#542788'),
    upset_query(set = 'ELF3', fill = '#8073ac'),
    upset_query(set = 'ELF4', fill = '#b2abd2')
  ),
  width_ratio = 0.1
)

#Compute intersection sizes (global)
count_all <- binary_df2 %>%
  dplyr::group_by(.intersection) %>%
  dplyr::count(.intersection) %>%
  arrange(desc(n))

top16 <- count_all$.intersection[1:16]

stack_order <- top16

tf_map <- c(
  "1" = "LHY",
  "2" = "CCA1",
  "3" = "TOC1",
  "4" = "PRR5",
  "5" = "PRR7",
  "6" = "LUX",
  "7" = "ELF3",
  "8" = "ELF4"
)

convert_intersection <- function(x) {
  parts <- strsplit(x, "-", fixed = TRUE)[[1]]
  paste(tf_map[parts], collapse = " & ")
}

top16_named <- sapply(top16, convert_intersection)

# Build amplitude stacked-bar annotations
amp_colors <- c(
  "gain high"   = "#31a354",
  "gain medium" = "#74c476",
  "other"       = "#cccccc",
  "lose medium" = "#fb6a4a",
  "lose high"   = "#de2d26"
)

# ---------------------------
# 1. Shared legend
# ---------------------------

legend_levels <- names(amp_colors)

legend_dummy <- ggplot(
  data.frame(
    x = 1,
    y = 1,
    fill = factor(legend_levels, levels = legend_levels)
  ),
  aes(x, y, fill = fill)
) +
  geom_col() +
  scale_fill_manual(values = amp_colors) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.background = element_rect(color = "black"),
    legend.box.spacing = unit(2.5, "pt")
  )

shared_legend <- cowplot::get_legend(legend_dummy)

# ---------------------------
# 2. d1d5 stacked bar (no x-axis labels)
# ---------------------------

stack_d1d5_plot <- stack_d1d5 +
  labs(title = "Amplitude: d1–d5") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

# ---------------------------
# 3. d1d2 stacked bar (with x-axis labels)
# ---------------------------

stack_d1d2_plot <- stack_d1d2 +
  labs(title = "Amplitude: d1–d2") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# ---------------------------
# 4. Panel A assembly
# ---------------------------

panel_A <- (
  wrap_elements(full = shared_legend) /
    plot_spacer() /
    stack_d1d5_plot /
    stack_d1d2_plot
) +
  plot_layout(heights = c(0.15, 0.05, 1, 1))

# ---------------------------
# 4. UpSet assembly (B)
# ---------------------------
up_panel <- wrap_elements(full = up)

# ---------------------------
# 5. Final figure: Panel A + Panel B (UpSet)
# ---------------------------

final_figure <- 
  wrap_elements(panel_A) / up_panel +
  plot_layout(heights = c(2.2, 3)) +
  plot_annotation(tag_levels = "A")

final_figure


#UpSet plot DTU-amp----
# Build CCA1‑Kamioka and CCA1‑Nagel subsets
kamioka_subset_DTU <- DTU_clock %>%
  filter(clock == 'CCA1-Kamioka') %>%
  mutate(clock = 'CCA1')

nagel_subset_DTU <- DTU_clock %>%
  filter(clock == 'CCA1-Nagel') %>%
  mutate(clock = 'CCA1')

# Intersection = true CCA1 DTU gene set
cca1_intersection_DTU <- inner_join(
  nagel_subset_DTU %>% select(gene_ID),
  kamioka_subset_DTU %>% select(gene_ID),
  by = "gene_ID",
  relationship = "many-to-many"   # optional, silences warning
) %>%
  distinct(gene_ID) %>%            # <-- THIS FIXES THE ISSUE
  mutate(clock = "CCA1")


# Build all DTU gene sets
LHY_DTU  <- DTU_clock %>% filter(clock == 'LHY')  %>% pull(gene_ID)
CCA1_DTU <- cca1_intersection_DTU %>% pull(gene_ID)
TOC1_DTU <- DTU_clock %>% filter(clock == 'TOC1') %>% pull(gene_ID)
PRR5_DTU <- DTU_clock %>% filter(clock == 'PRR5') %>% pull(gene_ID)
PRR7_DTU <- DTU_clock %>% filter(clock == 'PRR7') %>% pull(gene_ID)
LUX_DTU  <- DTU_clock %>% filter(clock == 'LUX')  %>% pull(gene_ID)
ELF3_DTU <- DTU_clock %>% filter(clock == 'ELF3') %>% pull(gene_ID)
ELF4_DTU <- DTU_clock %>% filter(clock == 'ELF4') %>% pull(gene_ID)

# Build the DTU gene‑set list
myGeneSets_DTU <- list(
  LHY  = LHY_DTU,
  CCA1 = CCA1_DTU,
  TOC1 = TOC1_DTU,
  PRR5 = PRR5_DTU,
  PRR7 = PRR7_DTU,
  LUX  = LUX_DTU,
  ELF3 = ELF3_DTU,
  ELF4 = ELF4_DTU
)

# Convert to long format
long_DTU <- enframe(myGeneSets_DTU, name = "clock", value = "gene_ID") %>%
  unnest(gene_ID) %>%
  distinct(clock, gene_ID)   # <-- critical fix

# Build binary matrix
binary_DTU <- long_DTU %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = clock,
    values_from = value,
    values_fill = 0
  )

# Build the DTU UpSet plot
# Build the UpSet plot FIRST (this defines the universe)
# due to bug with min_size move on to make intersection izes myself - see below
# use this to see the complete upset plot
# up_DTU <- ComplexUpset::upset(
#   binary_DTU,
#   intersect = c("LHY","CCA1","TOC1","PRR5","PRR7","LUX","ELF3","ELF4"),
#   #min_size = 5,
#   base_annotations = list(
#     "Intersection size" = intersection_size()
#   ),
#   set_sizes=(
#     upset_set_size()
#     + theme(axis.text.x=element_text(angle=90))
#   ),
#   queries = list(
#     upset_query(set = 'LHY',  fill = '#5aae61'),
#     upset_query(set = 'CCA1', fill = '#006d2c'),
#     upset_query(set = 'TOC1', fill = '#4d4d4d'),
#     upset_query(set = 'PRR5', fill = '#fdae61'),
#     upset_query(set = 'PRR7', fill = '#f46d43'),
#     upset_query(set = 'LUX',  fill = '#542788'),
#     upset_query(set = 'ELF3', fill = '#8073ac'),
#     upset_query(set = 'ELF4', fill = '#b2abd2')
#   ),
#   width_ratio = 0.1
# )

# Step 1 — compute intersection sizes myself
set_names <- c("LHY","CCA1","TOC1","PRR5","PRR7","LUX","ELF3","ELF4")

binary_DTU2 <- binary_DTU %>%
  rowwise() %>%
  mutate(
    .intersection = {
      present <- c_across(all_of(set_names)) == 1
      paste(set_names[present], collapse = "-")
    }
  ) %>%
  ungroup()

count_all_DTU <- binary_DTU2 %>%
  dplyr::group_by(.intersection) %>%
  dplyr::count(.intersection) %>%
  dplyr::arrange(desc(n))

top14 <- count_all_DTU$.intersection[count_all_DTU$n >= 2]

# Step 2 - filter the binary matrix to only those intersections
binary_DTU_filtered <- binary_DTU2 %>%
  filter(.intersection %in% top14)

# Step 3 - draw the UpSet plot without using min_size
up_DTU <- ComplexUpset::upset(
  binary_DTU_filtered,
  intersect = set_names,
  base_annotations = list(
    "Intersection size" = intersection_size()
  ),
  set_sizes = upset_set_size() +
    theme(axis.text.x = element_text(angle = 90)),
  queries = list(
    upset_query(set = 'LHY',  fill = '#5aae61'),
    upset_query(set = 'CCA1', fill = '#006d2c'),
    upset_query(set = 'TOC1', fill = '#4d4d4d'),
    upset_query(set = 'PRR5', fill = '#fdae61'),
    upset_query(set = 'PRR7', fill = '#f46d43'),
    upset_query(set = 'LUX',  fill = '#542788'),
    upset_query(set = 'ELF3', fill = '#8073ac'),
    upset_query(set = 'ELF4', fill = '#b2abd2')
  ),
  width_ratio = 0.1
)

# Build an amplitude metadata table
amp_meta_DTU <- DTU_d1d2_clock %>%
  select(gene_ID, amp_category) %>%
  distinct() %>%                     # <-- remove transcript duplicates
  rename(amp_d1d2 = amp_category) %>%
  full_join(
    DTU_d1d5_clock %>%
      select(gene_ID, amp_category) %>%
      distinct() %>%                 # <-- remove transcript duplicates
      rename(amp_d1d5 = amp_category),
    by = "gene_ID",
    relationship = "many-to-many"    # optional, silences warning
  ) %>%
  distinct(gene_ID, .keep_all = TRUE)  # <-- ensures one row per gene

# Build the DTU cluster‑level amplitude table
DTU_unique <- DTU_clock %>%
  distinct(gene_ID, cluster)

DTU_amp <- DTU_unique %>%
  left_join(amp_meta_DTU, by = "gene_ID")

# Order DTU clusters
cluster_order <- DTU_amp %>%
  dplyr::count(cluster) %>%
  dplyr::arrange(desc(n)) %>%
  pull(cluster)

DTU_amp_filtered <- DTU_amp %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

amp_levels <- c("gain high", "gain medium", "other", "lose medium", "lose high")

DTU_amp_filtered <- DTU_amp_filtered %>%
  mutate(
    amp_d1d5 = factor(amp_d1d5, levels = amp_levels),
    amp_d1d2 = factor(amp_d1d2, levels = amp_levels)
  )

# Build the DTU stacked bars
# DTU d1d5 stacked bar (top bar)
stack_DTU_d1d5 <- DTU_amp_filtered %>%
  dplyr::count(cluster, amp_d1d5) %>%
  ggplot(aes(x = cluster, y = n, fill = amp_d1d5)) +
  geom_col(position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = amp_colors,
                    breaks = amp_levels,
                    drop   = FALSE) +
  labs(
    x = "",
    y = "Proportion",
    title = "DTU Amplitude: d1–d5"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# DTU d1d2 stacked bar (bottom bar)
stack_DTU_d1d2 <- DTU_amp_filtered %>%
  dplyr::count(cluster, amp_d1d2) %>%
  ggplot(aes(x = cluster, y = n, fill = amp_d1d2)) +
  geom_col(position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = amp_colors,
                    breaks = amp_levels,
                    drop   = FALSE) +
  labs(
    x = "",
    y = "Proportion",
    title = "DTU Amplitude: d1–d2"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# Build DTU Panel A (stacked bars + legend)
# use the same shared legend in the DE figure for the DTU figure
# because the amplitude categories are identical

#shared_legend

panel_DTU_A <- (
  wrap_elements(full = shared_legend) /
    plot_spacer() /
    (stack_DTU_d1d5 + theme(plot.margin = margin())) /
    (stack_DTU_d1d2 + theme(plot.margin = margin()))
) +
  plot_layout(heights = c(0.15, 0.05, 1, 1))

panel_DTU_B <- wrap_elements(full = up_DTU)

final_DTU_figure <- 
  wrap_elements(panel_DTU_A) / panel_DTU_B +
  plot_layout(heights = c(2.2, 3)) +
  plot_annotation(tag_levels = "A")

# Rhythmicity stacked bar plots

# we have this from previous:
rhy_cols <- c(
  "nr-r"   = "#66c2a5",
  "r-nr"   = "#fc8d62",
  "r-r"    = "#006400",
  "nr-nr"  = "#bdbdbd"
)

# fixed order:
rhy_levels <- c("nr-r", "r-nr", "r-r", "nr-nr")

# Build the shared rhythmicity legend
legend_dummy_rhy <- ggplot(
  data.frame(
    x = 1,
    y = 1,
    fill = factor(rhy_levels, levels = rhy_levels)
  ),
  aes(x, y, fill = fill)
) +
  geom_col() +
  scale_fill_manual(values = rhy_cols, breaks = rhy_levels) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.box.background = element_rect(color = "black"),
    legend.box.spacing = unit(2.5, "pt")
  )

shared_rhy_legend <- cowplot::get_legend(legend_dummy_rhy)

# Build DE rhythmicity stacked bars
# Build an amplitude metadata table
rhy_meta_DE_d1d2 <- DE_d1d2_clock_rhy %>%
  select(gene_ID, rhy_category) 

rhy_meta_DE_d1d2 <- rhy_meta_DE_d1d2 %>%   
  mutate(
    rhy_category = factor(rhy_category, levels = rhy_levels)
  )

rhy_meta_DE_d1d5 <- DE_d1d5_clock_rhy %>%
  select(gene_ID, rhy_category) 

rhy_meta_DE_d1d5 <- rhy_meta_DE_d1d5 %>%   
  mutate(
    rhy_category = factor(rhy_category, levels = rhy_levels)
  )

# d1d5 stacked bar plot
stack_DE_rhy_d1d5 <- rhy_meta_DE_d1d5 %>%
  count(cluster, rhy_d1d5) %>%
  ggplot(aes(x = cluster, y = n, fill = rhy_d1d5)) +
  geom_col(position = "fill", width = 0.8) +
  scale_fill_manual(values = rhy_cols, breaks = rhy_levels, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Rhythmicity: d1–d5", x = "", y = "Proportion") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold")
  )













