source("scripts/00_config.R")
library(tidyverse); library(pheatmap); library(patchwork)
library(ggbeeswarm); library(RColorBrewer)

message("\n=== MODULE 5: Strain classification and clustering ===\n")

OUT <- file.path(RESULTS, "module05_classification")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

summary <- read_csv(file.path(DATA_PROC,"summary_metrics.csv"), show_col_types=FALSE)
meta    <- read_csv(FILE_METADATA, show_col_types=FALSE)
gmm     <- readRDS(file.path(RESULTS,"module02_exploratory","gmm_thresholds.rds"))

summary <- summary %>%
  mutate(across(c(max_delta_gamma,max_delta_gamma_OD,max_OD,
                  AUC_delta_gamma,AUC_delta_gamma_OD),
                ~if_else(is.infinite(.),NA_real_,.)),
         temperature = factor(temperature, levels=TEMP_LEVELS),
         salinity_m  = factor(salinity_m,  levels=c(0,0.5,1.0,1.5)))

message("GMM thresholds: Low/Mod=",round(gmm$thresh_low_moderate,1),
        " | Mod/High=",round(gmm$thresh_moderate_high,1))

# -----------------------------------------------------------------------------
# 1. Per-strain classification summary
# -----------------------------------------------------------------------------
message("\n--- Per-strain classification ---")

strain_class <- summary %>%
  filter(!is.na(max_delta_gamma_OD)) %>%
  group_by(strain) %>%
  summarise(
    max_ever         = max(max_delta_gamma_OD, na.rm=TRUE),
    mean_all         = mean(max_delta_gamma_OD, na.rm=TRUE),
    max_at_4C        = max(max_delta_gamma_OD[temperature=="4C"],  na.rm=TRUE),
    max_at_30C       = max(max_delta_gamma_OD[temperature=="30C"], na.rm=TRUE),
    max_AUC_4C       = max(AUC_delta_gamma_OD[temperature=="4C"],  na.rm=TRUE),
    max_AUC_30C      = max(AUC_delta_gamma_OD[temperature=="30C"], na.rm=TRUE),
    n_conditions     = n(),
    .groups="drop"
  ) %>%
  mutate(across(c(max_ever,mean_all,max_at_4C,max_at_30C,
                  max_AUC_4C,max_AUC_30C),
                ~if_else(is.infinite(.)|is.nan(.),NA_real_,.))) %>%
  mutate(
    producer_class = case_when(
      is.na(max_ever) | max_ever <= 0           ~ "Non-producer",
      max_ever < gmm$thresh_low_moderate         ~ "Low",
      max_ever < gmm$thresh_moderate_high        ~ "Moderate",
      TRUE                                       ~ "High"
    ),
    producer_class = factor(producer_class,
                            levels=c("High","Moderate","Low","Non-producer")),
    CSI = max_AUC_4C / max_AUC_30C,
    CSI_class = case_when(
      is.na(CSI)   ~ NA_character_,
      CSI > 1.5    ~ "Cold-stimulated",
      CSI < 0.7    ~ "Warm-preferred",
      TRUE         ~ "Neutral"
    )
  ) %>%
  left_join(meta %>% select(strain_phys,species,environment,
                             wgs_code,has_lipopeptide_bgc,bgc_type_lipopeptide),
            by=c("strain"="strain_phys"))

message("Producer classification:")
print(strain_class %>% count(producer_class))
message("\nCSI classification:")
print(strain_class %>% count(CSI_class))
message("\nFull strain classification table:")
print(strain_class %>%
        select(strain,producer_class,CSI_class,max_ever,
               max_at_4C,max_at_30C,species) %>%
        arrange(producer_class,desc(max_ever)) %>%
        mutate(across(where(is.numeric),~round(.,1))))

save_table(strain_class,"module05_strain_classification",dir=TABLES_MAIN)

# -----------------------------------------------------------------------------
# 2. Hierarchical clustering matrix
# -----------------------------------------------------------------------------
message("\n--- Hierarchical clustering ---")

clust_data <- summary %>%
  filter(!is.na(temperature), !is.na(max_delta_gamma_OD)) %>%
  mutate(condition=paste0(temperature,"_sal",as.character(salinity_m))) %>%
  group_by(strain,condition) %>%
  summarise(val=mean(max_delta_gamma_OD,na.rm=TRUE),.groups="drop") %>%
  pivot_wider(names_from=condition, values_from=val, values_fill=0) %>%
  column_to_rownames("strain")

message("Clustering matrix: ",nrow(clust_data)," strains x ",
        ncol(clust_data)," conditions")

# Ward D2 clustering
dist_mat  <- dist(scale(clust_data), method="euclidean")
hclust_obj <- hclust(dist_mat, method="ward.D2")
clusters   <- cutree(hclust_obj, k=3)
cluster_df <- data.frame(
  strain  = names(clusters),
  cluster = paste0("Cluster_",clusters)
)

message("Cluster assignments:")
print(cluster_df %>% left_join(strain_class %>%
        select(strain,producer_class,CSI_class),by="strain") %>%
        arrange(cluster))

save_table(cluster_df,"module05_clusters",dir=TABLES_MAIN)

# -----------------------------------------------------------------------------
# 3. OD vs production decoupling (formal Spearman test)
# -----------------------------------------------------------------------------
message("\n--- OD vs production decoupling ---")

decoupling <- summary %>%
  filter(!is.na(max_OD), !is.na(max_delta_gamma_OD),
         salinity_m==0) %>%
  group_by(temperature) %>%
  summarise(
    rho     = cor(max_OD, max_delta_gamma_OD, method="spearman"),
    n       = n(),
    p_value = cor.test(max_OD, max_delta_gamma_OD,
                       method="spearman")$p.value,
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric),~round(.,4)))

message("Spearman rho (OD vs production) per temperature:")
print(decoupling)
save_table(decoupling,"module05_decoupling_spearman",dir=TABLES_MAIN)

# -----------------------------------------------------------------------------
# 4. Figures
# -----------------------------------------------------------------------------
message("\nBuilding figures...")

# Fig A: Producer classification bar chart
p_producer <- strain_class %>%
  count(producer_class, environment) %>%
  filter(!is.na(producer_class)) %>%
  ggplot(aes(x=producer_class, y=n, fill=environment)) +
  geom_col(position="stack", width=0.6) +
  scale_fill_manual(values=COL_ENV, name="Origin",
                    na.value="grey70") +
  labs(x="Producer class", y="Number of strains",
       title="Biosurfactant producer classification",
       subtitle=paste0("GMM-based thresholds: Low/Mod=",
                       round(gmm$thresh_low_moderate,1),
                       " | Mod/High=",
                       round(gmm$thresh_moderate_high,1))) +
  THEME_PUB

# Fig B: CSI distribution plot
p_csi <- strain_class %>%
  filter(!is.na(CSI), !is.infinite(CSI)) %>%
  mutate(strain=fct_reorder(strain,CSI)) %>%
  ggplot(aes(x=CSI, y=strain, fill=CSI_class)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=c(0.7,1.5), linetype="dashed",
             colour="grey40", linewidth=0.5) +
  scale_fill_manual(values=COL_CSI, name="CSI class",
                    na.value="grey70") +
  scale_x_continuous(labels=scales::comma) +
  labs(x="Cold Stimulation Index (AUC 4C / AUC 30C)",
       y="Strain",
       title="Cold Stimulation Index per strain",
       subtitle="Dashed lines: thresholds at 0.7 (Warm) and 1.5 (Cold)") +
  THEME_PUB

# Fig C: Heatmap with clustering
ann_row <- strain_class %>%
  select(strain, producer_class, CSI_class, environment) %>%
  filter(strain %in% rownames(clust_data)) %>%
  distinct(strain, .keep_all=TRUE) %>%
  left_join(cluster_df, by="strain") %>%
  column_to_rownames("strain")

ann_colors <- list(
  producer_class = c(COL_PRODUCER,
                     "Non-producer"="grey90"),
  CSI_class      = c(COL_CSI, "Neutral"="#888780"),
  environment    = COL_ENV,
  cluster        = c(Cluster_1="#E8D5B7",
                     Cluster_2="#B7D5E8",
                     Cluster_3="#D5B7E8")
)

hmat <- as.matrix(clust_data)
hmat[is.na(hmat)] <- 0

out_heat <- file.path(OUT,"module05_clustering_heatmap.pdf")
pdf(out_heat, width=FIG_W_DOUBLE, height=FIG_H_MAX*0.7)
pheatmap(
  hmat,
  annotation_row   = ann_row[rownames(hmat),,drop=FALSE],
  annotation_colors = ann_colors,
  color            = colorRampPalette(
    c("white","#9FE1CB","#1D9E75","#085041"))(60),
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  fontsize         = 8,
  angle_col        = 45,
  main             = "Biosurfactant production across all conditions",
  na_col           = "grey95",
  border_color     = "white",
  scale            = "none"
)
dev.off()
message("  Heatmap saved: ", out_heat)

# Fig D: Max production at 4C vs 30C scatter
p_4v30 <- strain_class %>%
  filter(!is.na(max_at_4C), !is.na(max_at_30C),
         !is.infinite(max_at_4C), !is.infinite(max_at_30C)) %>%
  ggplot(aes(x=max_at_30C, y=max_at_4C,
             colour=CSI_class, label=strain)) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              colour="grey60", linewidth=0.6) +
  geom_point(size=3, alpha=0.85) +
  ggrepel::geom_text_repel(size=2.5, max.overlaps=12,
                            segment.colour="grey70",
                            segment.size=0.3) +
  scale_colour_manual(values=COL_CSI, name="CSI class",
                      na.value="grey70") +
  scale_x_continuous(labels=scales::comma) +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Max(dg/OD) at 30C",
       y="Max(dg/OD) at 4C",
       title="Cold vs warm biosurfactant efficiency",
       subtitle="Above diagonal = cold-stimulated | Below = warm-preferred") +
  THEME_PUB

if (!requireNamespace("ggrepel",quietly=TRUE))
  install.packages("ggrepel",repos="https://cloud.r-project.org")
library(ggrepel)

p_4v30 <- strain_class %>%
  filter(!is.na(max_at_4C),!is.na(max_at_30C),
         !is.infinite(max_at_4C),!is.infinite(max_at_30C)) %>%
  ggplot(aes(x=max_at_30C, y=max_at_4C,
             colour=CSI_class, label=strain)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",
              colour="grey60",linewidth=0.6) +
  geom_point(size=3,alpha=0.85) +
  geom_text_repel(size=2.5,max.overlaps=12,
                  segment.colour="grey70",segment.size=0.3) +
  scale_colour_manual(values=COL_CSI,name="CSI class",
                      na.value="grey70") +
  scale_x_continuous(labels=scales::comma) +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Max(dg/OD) at 30C", y="Max(dg/OD) at 4C",
       title="Cold vs warm biosurfactant efficiency",
       subtitle="Above diagonal = cold-stimulated") +
  THEME_PUB

# Combined figure
p_main <- (p_producer | p_csi) / (p_4v30 | plot_spacer()) +
  plot_annotation(title="Module 5: Strain classification",
                  tag_levels="A")

save_fig("module05_classification_main", p_main,
         w=FIG_W_DOUBLE, h=8, dir=OUT)
save_fig("module05_producer_classes",    p_producer, w=FIG_W_SINGLE, h=4, dir=OUT)
save_fig("module05_csi_ranking",         p_csi,      w=FIG_W_SINGLE, h=5, dir=OUT)
save_fig("module05_cold_vs_warm",        p_4v30,     w=FIG_W_SINGLE, h=4, dir=OUT)

# -----------------------------------------------------------------------------
# 5. Summary
# -----------------------------------------------------------------------------
message("\n=== MODULE 5 SUMMARY ===")
message("Total strains classified: ", nrow(strain_class))
message("\nProducer classification:")
print(strain_class %>% count(producer_class))
message("\nCSI classification:")
print(strain_class %>% count(CSI_class))
message("\nClusters (k=3, Ward D2):")
print(cluster_df %>% count(cluster))
message("\nDecoupling (OD vs production, Spearman rho):")
print(decoupling %>% select(temperature,rho,p_value))
message("\n=== MODULE 5 COMPLETE ===\n")
