source("scripts/00_config.R")
library(tidyverse); library(coin); library(patchwork); library(ggbeeswarm)

message("\n=== MODULE 6: Genomic integration ===\n")

OUT <- file.path(RESULTS, "module06_genomic")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
summary  <- read_csv(file.path(DATA_PROC,"summary_metrics.csv"), show_col_types=FALSE)
meta     <- read_csv(FILE_METADATA, show_col_types=FALSE)
strain_class <- read_csv(file.path(TABLES_MAIN,"module05_strain_classification.csv"),
                         show_col_types=FALSE)

# Load WGS inputs
bgc_matrix   <- read_tsv(FILE_BGC_MATRIX,   show_col_types=FALSE)
stress_genes <- read_tsv(FILE_STRESS_GENES, show_col_types=FALSE)
cog_matrix   <- read_tsv(FILE_COG_MATRIX,   show_col_types=FALSE)
wgs_meta     <- read_tsv(FILE_WGS_META,     show_col_types=FALSE)

message("BGC matrix: ",    nrow(bgc_matrix),   " rows x ", ncol(bgc_matrix),   " cols")
message("Stress genes: ",  nrow(stress_genes), " rows x ", ncol(stress_genes), " cols")
message("COG matrix: ",    nrow(cog_matrix),   " rows x ", ncol(cog_matrix),   " cols")
message("WGS metadata: ",  nrow(wgs_meta),     " rows x ", ncol(wgs_meta),     " cols")

message("\nBGC matrix columns: ",  paste(names(bgc_matrix),  collapse=" | "))
message("Stress genes columns: ", paste(names(stress_genes),collapse=" | "))
message("WGS meta columns: ",     paste(names(wgs_meta),    collapse=" | "))

# -----------------------------------------------------------------------------
# 2. Build integrated dataset
# -----------------------------------------------------------------------------
message("\n--- Building integrated genotype-phenotype dataset ---")

# Physiology summary per strain at 4C (primary condition of interest)
phys_4C <- summary %>%
  filter(!is.na(max_delta_gamma_OD), !is.na(temperature)) %>%
  filter(temperature=="4C") %>%
  group_by(strain) %>%
  summarise(
    max_dg_OD_4C  = max(max_delta_gamma_OD, na.rm=TRUE),
    mean_dg_OD_4C = mean(max_delta_gamma_OD, na.rm=TRUE),
    AUC_4C        = max(AUC_delta_gamma_OD,  na.rm=TRUE),
    max_OD_4C     = max(max_OD,              na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric), ~if_else(is.infinite(.)|is.nan(.),NA_real_,.)))

# Join with metadata
gp_data <- phys_4C %>%
  left_join(meta %>% select(strain_phys, wgs_code, species, environment,
                             has_lipopeptide_bgc, bgc_type_lipopeptide,
                             has_wgs, CSI_class),
            by=c("strain"="strain_phys")) %>%
  left_join(strain_class %>% select(strain, producer_class, max_ever),
            by="strain") %>%
  filter(has_wgs == TRUE)

message("Strains with both physiology and WGS: ", nrow(gp_data))
message("Strains: ", paste(gp_data$strain, collapse=", "))

# -----------------------------------------------------------------------------
# 3. H3: Lipopeptide BGC presence vs Max(dg/OD) at 4C
# -----------------------------------------------------------------------------
message("\n--- H3: Lipopeptide BGC vs production at 4C ---")

bgc_pos <- gp_data %>% filter(has_lipopeptide_bgc==TRUE)
bgc_neg <- gp_data %>% filter(has_lipopeptide_bgc==FALSE | is.na(has_lipopeptide_bgc))

message("BGC-positive strains: ", nrow(bgc_pos), " — ",
        paste(bgc_pos$strain, collapse=", "))
message("BGC-negative strains: ", nrow(bgc_neg), " — ",
        paste(bgc_neg$strain, collapse=", "))

message("\nBGC-positive: median Max(dg/OD) at 4C = ",
        round(median(bgc_pos$max_dg_OD_4C, na.rm=TRUE),1))
message("BGC-negative: median Max(dg/OD) at 4C = ",
        round(median(bgc_neg$max_dg_OD_4C, na.rm=TRUE),1))

# Wilcoxon rank-sum test
wt <- wilcox.test(max_dg_OD_4C ~ has_lipopeptide_bgc,
                  data=gp_data %>% filter(!is.na(has_lipopeptide_bgc)),
                  exact=FALSE)
message("Wilcoxon test: W=",round(wt$statistic,1)," p=",round(wt$p.value,4))

# Effect size: rank-biserial correlation
n1 <- nrow(bgc_pos %>% filter(!is.na(max_dg_OD_4C)))
n2 <- nrow(bgc_neg %>% filter(!is.na(max_dg_OD_4C)))
r_rb <- 1 - (2*wt$statistic)/(n1*n2)
message("Rank-biserial correlation (effect size): r=",round(r_rb,3))

bgc_test_result <- data.frame(
  test="Wilcoxon rank-sum",
  W=wt$statistic, p_value=wt$p.value,
  n_bgc_pos=n1, n_bgc_neg=n2,
  median_bgc_pos=median(bgc_pos$max_dg_OD_4C,na.rm=TRUE),
  median_bgc_neg=median(bgc_neg$max_dg_OD_4C,na.rm=TRUE),
  rank_biserial_r=r_rb
) %>% mutate(across(where(is.numeric),~round(.,4)))
save_table(bgc_test_result,"module06_bgc_wilcoxon",dir=TABLES_MAIN)

# -----------------------------------------------------------------------------
# 4. BGC count vs production (Spearman)
# -----------------------------------------------------------------------------
message("\n--- BGC count vs production ---")

# Get BGC counts from WGS metadata
bgc_count_col <- names(wgs_meta)[str_detect(names(wgs_meta),
                                  regex("bgc.*count|n_bgc|total.*bgc",
                                        ignore_case=TRUE))]
message("BGC count column: ", paste(bgc_count_col, collapse=", "))

# Also get strain ID column
strain_col_wgs <- names(wgs_meta)[str_detect(names(wgs_meta),
                                   regex("strain|isolate|sample",
                                         ignore_case=TRUE))][1]
message("Strain col in WGS meta: ", strain_col_wgs)
message("First few rows of WGS meta:")
print(head(wgs_meta[,1:min(6,ncol(wgs_meta))]))

# -----------------------------------------------------------------------------
# 5. Stress gene analysis
# -----------------------------------------------------------------------------
message("\n--- Stress gene analysis ---")
message("Stress genes columns:")
print(names(stress_genes))
message("First few rows:")
print(head(stress_genes))

# -----------------------------------------------------------------------------
# 6. Surfactin vs Lichenysin producers
# -----------------------------------------------------------------------------
message("\n--- Surfactin vs Lichenysin producers ---")

lipopeptide_comp <- gp_data %>%
  filter(!is.na(bgc_type_lipopeptide)) %>%
  group_by(bgc_type_lipopeptide) %>%
  summarise(
    n             = n(),
    median_dg_OD  = median(max_dg_OD_4C, na.rm=TRUE),
    mean_dg_OD    = mean(max_dg_OD_4C,   na.rm=TRUE),
    median_AUC    = median(AUC_4C,        na.rm=TRUE),
    strains       = paste(strain, collapse=", "),
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric),~round(.,2)))

message("Surfactin vs Lichenysin comparison:")
print(lipopeptide_comp)
save_table(lipopeptide_comp,"module06_lipopeptide_comparison",dir=TABLES_MAIN)

# -----------------------------------------------------------------------------
# 7. H2E-10C paradox analysis
# -----------------------------------------------------------------------------
message("\n--- H2E-10C paradox: cold-stimulated without lipopeptide BGC ---")
h2e_data <- gp_data %>% filter(strain=="H2E-10C")
message("H2E-10C: WGS=",h2e_data$wgs_code,
        " | Lipopeptide BGC=",h2e_data$has_lipopeptide_bgc,
        " | Max(dg/OD) at 4C=",round(h2e_data$max_dg_OD_4C,1),
        " | CSI=",h2e_data$CSI_class)
message("Species: P. frigoritolerans (cold-tolerant by name)")
message("This strain produces biosurfactant without canonical lipopeptide BGC")
message("Check antiSMASH output at: /data/habib/WGS/analysis/06_antismash/A6/")

# -----------------------------------------------------------------------------
# 8. Environment comparison: Qaidam vs Marine
# -----------------------------------------------------------------------------
message("\n--- Qaidam vs Marine production at 4C ---")

env_comp <- gp_data %>%
  filter(!is.na(environment), !is.na(max_dg_OD_4C)) %>%
  group_by(environment) %>%
  summarise(
    n            = n(),
    median_dg_OD = median(max_dg_OD_4C, na.rm=TRUE),
    mean_dg_OD   = mean(max_dg_OD_4C,   na.rm=TRUE),
    sd_dg_OD     = sd(max_dg_OD_4C,     na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric),~round(.,2)))
message("Environment comparison:")
print(env_comp)

env_wt <- wilcox.test(max_dg_OD_4C ~ environment,
                      data=gp_data %>% filter(!is.na(environment),
                                              !is.na(max_dg_OD_4C)),
                      exact=FALSE)
message("Wilcoxon Qaidam vs Marine: W=",round(env_wt$statistic,1),
        " p=",round(env_wt$p.value,4))

save_table(env_comp,"module06_environment_comparison",dir=TABLES_MAIN)

# -----------------------------------------------------------------------------
# 9. Figures
# -----------------------------------------------------------------------------
message("\nBuilding figures...")

# Fig A: BGC presence vs production at 4C
p_bgc <- gp_data %>%
  filter(!is.na(has_lipopeptide_bgc), !is.na(max_dg_OD_4C)) %>%
  mutate(bgc_label=if_else(has_lipopeptide_bgc,
                            "BGC-positive\n(Lichenysin/Surfactin)",
                            "BGC-negative")) %>%
  ggplot(aes(x=bgc_label, y=max_dg_OD_4C, colour=bgc_label)) +
  geom_boxplot(aes(fill=bgc_label), alpha=0.15, outlier.shape=NA,
               linewidth=0.5, width=0.45) +
  geom_beeswarm(size=2.5, cex=2.5, alpha=0.85) +
  ggrepel::geom_text_repel(aes(label=strain), size=2.2,
                            max.overlaps=10, segment.size=0.3,
                            segment.colour="grey70") +
  scale_colour_manual(values=c("BGC-positive\n(Lichenysin/Surfactin)"="#1D9E75",
                                "BGC-negative"="#D85A30"), guide="none") +
  scale_fill_manual(values=c("BGC-positive\n(Lichenysin/Surfactin)"="#1D9E75",
                              "BGC-negative"="#D85A30"), guide="none") +
  scale_y_continuous(labels=scales::comma) +
  annotate("text", x=1.5, y=max(gp_data$max_dg_OD_4C,na.rm=TRUE)*0.95,
           label=paste0("W=",round(wt$statistic,0),
                        "\np=",round(wt$p.value,3)),
           size=3, colour="grey30") +
  labs(x="Lipopeptide BGC status",
       y="Max(dg/OD) at 4C",
       title="H3: Lipopeptide BGC predicts cold biosurfactant production",
       subtitle=paste0("Wilcoxon W=",round(wt$statistic,1),
                       " p=",round(wt$p.value,4),
                       " | r=",round(r_rb,3))) +
  THEME_PUB

library(ggrepel)
p_bgc <- gp_data %>%
  filter(!is.na(has_lipopeptide_bgc), !is.na(max_dg_OD_4C)) %>%
  mutate(bgc_label=if_else(has_lipopeptide_bgc,
                            "BGC-positive","BGC-negative")) %>%
  ggplot(aes(x=bgc_label, y=max_dg_OD_4C, colour=bgc_label)) +
  geom_boxplot(aes(fill=bgc_label), alpha=0.15, outlier.shape=NA,
               linewidth=0.5, width=0.45) +
  geom_beeswarm(size=2.5, cex=2.5, alpha=0.85) +
  geom_text_repel(aes(label=strain), size=2.2, max.overlaps=10,
                  segment.size=0.3, segment.colour="grey70") +
  scale_colour_manual(values=c("BGC-positive"="#1D9E75",
                                "BGC-negative"="#D85A30"), guide="none") +
  scale_fill_manual(values=c("BGC-positive"="#1D9E75",
                              "BGC-negative"="#D85A30"), guide="none") +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Lipopeptide BGC status", y="Max(dg/OD) at 4C",
       title="H3: Lipopeptide BGC predicts cold biosurfactant production",
       subtitle=paste0("Wilcoxon W=",round(wt$statistic,1),
                       " p=",round(wt$p.value,4),
                       " r=",round(r_rb,3))) +
  THEME_PUB

# Fig B: Surfactin vs Lichenysin
p_lipopeptide <- gp_data %>%
  filter(!is.na(bgc_type_lipopeptide), !is.na(max_dg_OD_4C)) %>%
  ggplot(aes(x=bgc_type_lipopeptide, y=max_dg_OD_4C,
             colour=bgc_type_lipopeptide)) +
  geom_boxplot(aes(fill=bgc_type_lipopeptide), alpha=0.15,
               outlier.shape=NA, linewidth=0.5, width=0.45) +
  geom_beeswarm(size=2.5, cex=2.5, alpha=0.85) +
  geom_text_repel(aes(label=strain), size=2.2, max.overlaps=10,
                  segment.size=0.3, segment.colour="grey70") +
  scale_colour_manual(values=c("Surfactin"="#534AB7",
                                "Lichenysin"="#1D9E75"), guide="none") +
  scale_fill_manual(values=c("Surfactin"="#534AB7",
                              "Lichenysin"="#1D9E75"), guide="none") +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Lipopeptide type", y="Max(dg/OD) at 4C",
       title="Surfactin vs Lichenysin producers at 4C",
       subtitle="SW2C3 is the only Surfactin producer") +
  THEME_PUB

# Fig C: Qaidam vs Marine
p_env <- gp_data %>%
  filter(!is.na(environment), !is.na(max_dg_OD_4C)) %>%
  ggplot(aes(x=environment, y=max_dg_OD_4C, colour=environment)) +
  geom_boxplot(aes(fill=environment), alpha=0.15, outlier.shape=NA,
               linewidth=0.5, width=0.45) +
  geom_beeswarm(size=2.5, cex=2.5, alpha=0.85) +
  geom_text_repel(aes(label=strain), size=2.2, max.overlaps=10,
                  segment.size=0.3, segment.colour="grey70") +
  scale_colour_manual(values=COL_ENV, guide="none") +
  scale_fill_manual(values=COL_ENV,   guide="none") +
  scale_y_continuous(labels=scales::comma) +
  annotate("text", x=1.5,
           y=max(gp_data$max_dg_OD_4C,na.rm=TRUE)*0.95,
           label=paste0("p=",round(env_wt$p.value,3)),
           size=3, colour="grey30") +
  labs(x="Origin", y="Max(dg/OD) at 4C",
       title="Qaidam vs Marine production at 4C",
       subtitle=paste0("Wilcoxon p=",round(env_wt$p.value,4))) +
  THEME_PUB

# Fig D: CSI vs BGC presence
p_csi_bgc <- gp_data %>%
  filter(!is.na(has_lipopeptide_bgc), !is.na(CSI_class)) %>%
  mutate(bgc_label=if_else(has_lipopeptide_bgc,
                            "BGC-positive","BGC-negative")) %>%
  count(bgc_label, CSI_class) %>%
  ggplot(aes(x=bgc_label, y=n, fill=CSI_class)) +
  geom_col(position="fill", width=0.6) +
  scale_fill_manual(values=COL_CSI, name="CSI class") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Lipopeptide BGC status",
       y="Proportion of strains",
       title="CSI class by BGC status",
       subtitle="BGC-positive strains enriched in Cold-stimulated class") +
  THEME_PUB

p_main <- (p_bgc | p_lipopeptide) / (p_env | p_csi_bgc) +
  plot_annotation(title="Module 6: Genomic integration",
                  tag_levels="A")

save_fig("module06_genomic_main",     p_main,        w=FIG_W_DOUBLE, h=8, dir=OUT)
save_fig("module06_bgc_vs_production",p_bgc,         w=FIG_W_SINGLE, h=4, dir=OUT)
save_fig("module06_lipopeptide_type", p_lipopeptide, w=FIG_W_SINGLE, h=4, dir=OUT)
save_fig("module06_environment",      p_env,         w=FIG_W_SINGLE, h=4, dir=OUT)
save_fig("module06_csi_by_bgc",       p_csi_bgc,     w=FIG_W_SINGLE, h=4, dir=OUT)

# -----------------------------------------------------------------------------
# 10. Summary
# -----------------------------------------------------------------------------
message("\n=== MODULE 6 SUMMARY ===")
message("Strains with WGS + physiology: ", nrow(gp_data))
message("BGC-positive (lipopeptide):    ", sum(gp_data$has_lipopeptide_bgc, na.rm=TRUE))
message("BGC-negative:                  ", sum(!gp_data$has_lipopeptide_bgc, na.rm=TRUE))
message("Wilcoxon BGC vs production:    W=",round(wt$statistic,1),
        " p=",round(wt$p.value,4))
message("Effect size (rank-biserial r): ", round(r_rb,3))
message("Qaidam vs Marine Wilcoxon p:   ", round(env_wt$p.value,4))
message("H2E-10C: paradox strain noted (cold-stimulated, no lipopeptide BGC)")
message("\n=== MODULE 6 COMPLETE ===\n")
