source("scripts/00_config.R")
library(tidyverse); library(mclust); library(pheatmap); library(ggbeeswarm)
library(patchwork); library(RColorBrewer)

message("\n=== MODULE 2: Exploratory analysis ===\n")

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
master  <- read_csv(file.path(DATA_PROC, "master_dataset.csv"),  show_col_types=FALSE)
summary <- read_csv(file.path(DATA_PROC, "summary_metrics.csv"), show_col_types=FALSE)
meta    <- read_csv(FILE_METADATA, show_col_types=FALSE)

# Clean -Inf from max() on empty groups
summary <- summary %>%
  mutate(across(c(max_delta_gamma, max_delta_gamma_OD, max_OD,
                  AUC_delta_gamma, AUC_delta_gamma_OD),
                ~if_else(is.infinite(.), NA_real_, .)))

master <- master %>%
  mutate(temperature = factor(temperature, levels=TEMP_LEVELS),
         salinity_m  = factor(salinity_m,  levels=c(0, 0.5, 1.0, 1.5)))

summary <- summary %>%
  mutate(temperature = factor(temperature, levels=TEMP_LEVELS),
         salinity_m  = factor(salinity_m,  levels=c(0, 0.5, 1.0, 1.5)))

message("Data loaded: ", nrow(master), " readings | ",
        nrow(summary), " condition summaries")

# -----------------------------------------------------------------------------
# 2. Descriptive statistics table
# -----------------------------------------------------------------------------
message("\nComputing descriptive statistics...")

desc_stats <- summary %>%
  group_by(temperature, salinity_m) %>%
  summarise(
    n_strains        = sum(!is.na(max_delta_gamma_OD)),
    mean_max_dg_OD   = mean(max_delta_gamma_OD, na.rm=TRUE),
    sd_max_dg_OD     = sd(max_delta_gamma_OD,   na.rm=TRUE),
    median_max_dg_OD = median(max_delta_gamma_OD, na.rm=TRUE),
    min_max_dg_OD    = min(max_delta_gamma_OD,  na.rm=TRUE),
    max_max_dg_OD    = max(max_delta_gamma_OD,  na.rm=TRUE),
    mean_AUC         = mean(AUC_delta_gamma_OD, na.rm=TRUE),
    sd_AUC           = sd(AUC_delta_gamma_OD,   na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

save_table(desc_stats, "module02_descriptive_stats", dir=TABLES_MAIN)
message("  Descriptive stats saved")
print(desc_stats %>% filter(salinity_m==0) %>% select(temperature, n_strains,
      mean_max_dg_OD, sd_max_dg_OD, median_max_dg_OD))

# -----------------------------------------------------------------------------
# 3. Gaussian mixture model — producer classification
# -----------------------------------------------------------------------------
message("\nFitting Gaussian mixture model for producer classification...")

# Use Max(delta_gamma_OD) pooled across ALL conditions for classification
all_vals <- summary %>%
  filter(!is.na(max_delta_gamma_OD), max_delta_gamma_OD > 0) %>%
  pull(max_delta_gamma_OD)

message("  Values for GMM: n=", length(all_vals),
        " range=[", round(min(all_vals),1), ", ", round(max(all_vals),1), "]")

# Fit GMM with BIC selection (1-4 components)
gmm_fit <- Mclust(all_vals, G=1:4, verbose=FALSE)
message("  BIC-selected components: ", gmm_fit$G)
message("  Component means: ", paste(round(sort(gmm_fit$parameters$mean),1), collapse=", "))

# Get thresholds as intersection points between components
# For manuscript: use the component means to define thresholds
comp_means <- sort(gmm_fit$parameters$mean)
if (gmm_fit$G >= 2) {
  thresh_low  <- mean(comp_means[1:2])
  thresh_high <- if (gmm_fit$G >= 3) mean(comp_means[(gmm_fit$G-1):gmm_fit$G]) else thresh_low
} else {
  thresh_low  <- median(all_vals)
  thresh_high <- quantile(all_vals, 0.75)
}
message("  Thresholds: Low/Moderate=", round(thresh_low,1),
        " | Moderate/High=", round(thresh_high,1))

# Assign classification to summary
summary <- summary %>%
  mutate(
    producer_class = case_when(
      is.na(max_delta_gamma_OD) | max_delta_gamma_OD <= 0 ~ "Non-producer",
      max_delta_gamma_OD < thresh_low                      ~ "Low",
      max_delta_gamma_OD < thresh_high                     ~ "Moderate",
      TRUE                                                 ~ "High"
    ),
    producer_class = factor(producer_class,
                            levels=c("High","Moderate","Low","Non-producer"))
  )

prod_summary <- summary %>%
  filter(!is.na(producer_class)) %>%
  group_by(strain) %>%
  slice_max(order_by=as.integer(producer_class), n=1, with_ties=FALSE) %>%
  ungroup() %>%
  count(producer_class)
message("  Producer classification:")
print(prod_summary)

# Save GMM thresholds for reference in other modules
gmm_thresholds <- list(
  n_components = gmm_fit$G,
  component_means = comp_means,
  thresh_low_moderate = thresh_low,
  thresh_moderate_high = thresh_high
)
saveRDS(gmm_thresholds, file.path(RESULTS, "module02_exploratory", "gmm_thresholds.rds"))

# Update summary file with producer class
write_csv(summary, file.path(DATA_PROC, "summary_metrics.csv"))

# -----------------------------------------------------------------------------
# 4. Figure 1 — GMM distribution plot
# -----------------------------------------------------------------------------
message("\nBuilding figures...")

p_gmm <- ggplot(data.frame(x=all_vals), aes(x=x)) +
  geom_histogram(aes(y=after_stat(density)), bins=40,
                 fill="grey80", colour="grey60", linewidth=0.3) +
  geom_density(colour="grey40", linewidth=0.6) +
  geom_vline(xintercept=thresh_low,  linetype="dashed",
             colour="#EF9F27", linewidth=0.8) +
  geom_vline(xintercept=thresh_high, linetype="dashed",
             colour="#1D9E75", linewidth=0.8) +
  annotate("text", x=thresh_low*0.5,  y=Inf, vjust=1.5, size=3,
           label="Low", colour="#EF9F27") +
  annotate("text", x=(thresh_low+thresh_high)/2, y=Inf, vjust=1.5, size=3,
           label="Moderate", colour="grey40") +
  annotate("text", x=thresh_high*1.3, y=Inf, vjust=1.5, size=3,
           label="High", colour="#1D9E75") +
  labs(x=expression("Max("*Delta*gamma*"/OD) (mN/m per OD)"),
       y="Density",
       title="Producer classification via Gaussian mixture model",
       subtitle=paste0("BIC-selected components: ", gmm_fit$G,
                       " | n=", length(all_vals), " condition values")) +
  THEME_PUB

save_fig("module02_gmm_distribution", p_gmm,
         w=FIG_W_SINGLE, h=3.5,
         dir=file.path(RESULTS,"module02_exploratory"))

# -----------------------------------------------------------------------------
# 5. Figure 2 — Baseline ST homogeneity across experiments
# -----------------------------------------------------------------------------
baseline_data <- master %>%
  filter(day==1) %>%
  mutate(experiment=factor(experiment,
         levels=c("exp2","exp3","exp4","exp5","exp6",
                  "exp7","exp8","exp9","exp10","exp11")))

p_baseline <- ggplot(baseline_data, aes(x=experiment, y=st)) +
  geom_boxplot(fill="grey90", colour="grey40", outlier.shape=NA,
               linewidth=0.4, width=0.5) +
  geom_jitter(aes(colour=temperature), width=0.15, size=1.5, alpha=0.7) +
  geom_hline(yintercept=c(68,71), linetype="dotted",
             colour="red", linewidth=0.5) +
  scale_colour_manual(values=COL_TEMP, name="Temperature") +
  scale_y_continuous(limits=c(67,72), breaks=seq(67,72,1)) +
  labs(x="Experiment", y="Baseline ST (mN/m)",
       title="Baseline surface tension homogeneity",
       subtitle="Day 1 readings across all experiments (dotted lines = 68-71 mN/m expected range)") +
  THEME_PUB +
  theme(axis.text.x=element_text(angle=45, hjust=1))

save_fig("module02_baseline_homogeneity", p_baseline,
         w=FIG_W_DOUBLE, h=3.5,
         dir=file.path(RESULTS,"module02_exploratory"))

# -----------------------------------------------------------------------------
# 6. Figure 3 — Max(delta_gamma_OD) by temperature, all strains
# -----------------------------------------------------------------------------
temp_only_summary <- summary %>%
  filter(as.numeric(as.character(salinity_m))==0, !is.na(max_delta_gamma_OD), !is.na(temperature))

p_temp_overview <- ggplot(temp_only_summary,
                          aes(x=temperature, y=max_delta_gamma_OD,
                              colour=temperature)) +
  geom_boxplot(aes(fill=temperature), alpha=0.15, outlier.shape=NA,
               linewidth=0.5, width=0.5) +
  geom_beeswarm(size=2, cex=2.5, alpha=0.85) +
  scale_colour_manual(values=COL_TEMP, guide="none") +
  scale_fill_manual(values=COL_TEMP,   guide="none") +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Temperature", y=expression("Max("*Delta*gamma*"/OD)"),
       title="Biosurfactant efficiency by temperature",
       subtitle="All strains, temp-only experiments (salinity=0)") +
  THEME_PUB

save_fig("module02_temp_overview", p_temp_overview,
         w=FIG_W_SINGLE, h=3.5,
         dir=file.path(RESULTS,"module02_exploratory"))

# -----------------------------------------------------------------------------
# 7. Figure 4 — OD vs delta_gamma_OD scatter (growth vs production decoupling)
# -----------------------------------------------------------------------------
p_decoupling <- ggplot(
  temp_only_summary %>% filter(!is.na(max_OD), !is.na(max_delta_gamma_OD)),
  aes(x=max_OD, y=max_delta_gamma_OD,
      colour=temperature, label=strain)) +
  geom_point(size=2.5, alpha=0.85) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.6, alpha=0.15) +
  ggrepel::geom_text_repel(size=2.5, max.overlaps=8,
                            segment.colour="grey70") +
  scale_colour_manual(values=COL_TEMP, name="Temperature") +
  scale_y_continuous(labels=scales::comma) +
  facet_wrap(~temperature, nrow=1) +
  labs(x="Max OD (growth)", y=expression("Max("*Delta*gamma*"/OD)"),
       title="Growth vs biosurfactant production",
       subtitle="No correlation = production is decoupled from biomass") +
  THEME_PUB +
  theme(legend.position="none")

# Need ggrepel
if (!requireNamespace("ggrepel", quietly=TRUE)) {
  message("  Installing ggrepel...")
  install.packages("ggrepel", repos="https://cloud.r-project.org")
}
library(ggrepel)

p_decoupling <- ggplot(
  temp_only_summary %>% filter(!is.na(max_OD), !is.na(max_delta_gamma_OD)),
  aes(x=max_OD, y=max_delta_gamma_OD, colour=temperature, label=strain)) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.6, alpha=0.12,
              colour="grey50") +
  geom_point(size=2.5, alpha=0.85) +
  geom_text_repel(size=2.2, max.overlaps=10, segment.colour="grey70",
                  segment.size=0.3) +
  scale_colour_manual(values=COL_TEMP, name="Temperature") +
  scale_y_continuous(labels=scales::comma) +
  facet_wrap(~temperature, nrow=1) +
  labs(x="Max OD (growth proxy)",
       y=expression("Max("*Delta*gamma*"/OD) (biosurfactant efficiency)"),
       title="Production is decoupled from growth",
       subtitle="Spearman correlation per temperature shown") +
  THEME_PUB + theme(legend.position="none")

save_fig("module02_growth_vs_production", p_decoupling,
         w=FIG_W_DOUBLE, h=3.8,
         dir=file.path(RESULTS,"module02_exploratory"))

# Compute Spearman correlations
message("\nSpearman correlation: Max(OD) vs Max(delta_gamma_OD) per temperature:")
temp_only_summary %>%
  filter(!is.na(max_OD), !is.na(max_delta_gamma_OD)) %>%
  group_by(temperature) %>%
  summarise(
    rho   = cor(max_OD, max_delta_gamma_OD, method="spearman"),
    n     = n(),
    .groups="drop"
  ) %>%
  mutate(rho=round(rho,3)) %>%
  print()

# -----------------------------------------------------------------------------
# 8. Figure 5 — Heatmap: all strains × temperature (temp-only)
# -----------------------------------------------------------------------------
message("\nBuilding heatmap...")

heatmap_data <- temp_only_summary %>%
  filter(!is.na(temperature)) %>%
  group_by(strain, temperature) %>%
  summarise(val=mean(max_delta_gamma_OD, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=temperature, values_from=val) %>%
  column_to_rownames("strain")

# Order columns correctly
col_order <- intersect(TEMP_LEVELS, names(heatmap_data))
heatmap_data <- heatmap_data[, col_order, drop=FALSE]

# Annotation: CSI class and environment
ann_row <- meta %>%
  filter(strain_phys %in% rownames(heatmap_data)) %>%
  select(strain_phys, CSI_class, environment) %>%
  distinct(strain_phys, .keep_all=TRUE) %>%
  column_to_rownames("strain_phys")

ann_colors <- list(
  CSI_class   = COL_CSI,
  environment = COL_ENV
)

heatmap_matrix <- as.matrix(heatmap_data)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

out_heatmap <- file.path(RESULTS, "module02_exploratory",
                         "module02_heatmap_all_strains.pdf")
dir.create(dirname(out_heatmap), recursive=TRUE, showWarnings=FALSE)

pdf(out_heatmap, width=FIG_W_SINGLE*1.5, height=FIG_H_MAX*0.6)
pheatmap(heatmap_matrix,
         annotation_row = ann_row[rownames(heatmap_matrix),, drop=FALSE],
         annotation_colors = ann_colors,
         color = colorRampPalette(c("white","#9FE1CB","#1D9E75","#085041"))(50),
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 9,
         main = "Max(Δγ/OD) across temperatures — all strains",
         na_col = "grey95",
         border_color = "white")
dev.off()
message("  Heatmap saved: ", out_heatmap)

# -----------------------------------------------------------------------------
# 9. Summary report
# -----------------------------------------------------------------------------
message("\n=== MODULE 2 SUMMARY ===")
message("GMM components:       ", gmm_fit$G)
message("Threshold Low/Mod:    ", round(thresh_low,1))
message("Threshold Mod/High:   ", round(thresh_high,1))
message("Producer classes:")
print(prod_summary)
message("Figures saved to:     results/module02_exploratory/")
message("Tables saved to:      tables/main/")
message("\n=== MODULE 2 COMPLETE ===\n")
