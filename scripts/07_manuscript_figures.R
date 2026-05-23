source("scripts/00_config.R")
library(tidyverse); library(patchwork); library(ggbeeswarm)
library(pheatmap); library(RColorBrewer); library(ggrepel)
library(lme4); library(lmerTest); library(emmeans)

message("\n=== MODULE 7: Manuscript figures ===\n")

OUT_MAIN <- FIGS_MAIN
OUT_SUPP <- FIGS_SUPP
dir.create(OUT_MAIN, recursive=TRUE, showWarnings=FALSE)
dir.create(OUT_SUPP, recursive=TRUE, showWarnings=FALSE)

# Load all processed data
master       <- read_csv(file.path(DATA_PROC,"master_dataset.csv"),    show_col_types=FALSE)
summary_df   <- read_csv(file.path(DATA_PROC,"summary_metrics.csv"),   show_col_types=FALSE)
meta         <- read_csv(FILE_METADATA, show_col_types=FALSE)
strain_class <- read_csv(file.path(TABLES_MAIN,"module05_strain_classification.csv"), show_col_types=FALSE)
bliss        <- read_csv(file.path(TABLES_MAIN,"module04_bliss_synergy.csv"),         show_col_types=FALSE)
stress_spear <- read_csv(file.path(TABLES_MAIN,"module06_stress_gene_spearman.csv"),  show_col_types=FALSE)
strain_prof  <- read_csv(file.path(TABLES_MAIN,"module04_strain_profiles.csv"),       show_col_types=FALSE)

# Clean and factor
summary_df <- summary_df %>%
  mutate(across(c(max_delta_gamma,max_delta_gamma_OD,max_OD,
                  AUC_delta_gamma,AUC_delta_gamma_OD),
                ~if_else(is.infinite(.),NA_real_,.)),
         temperature = factor(temperature, levels=TEMP_LEVELS),
         salinity_m  = factor(salinity_m,  levels=c(0,0.5,1.0,1.5)))

master <- master %>%
  mutate(temperature = factor(temperature, levels=TEMP_LEVELS),
         salinity_m  = factor(salinity_m,  levels=c(0,0.5,1.0,1.5)))

message("Data loaded successfully")

# =============================================================================
# FIGURE 1: The Cold-Induction Phenomenon
# Time-series of key cold-stimulated strains at 4C vs 30C
# =============================================================================
message("\nBuilding Figure 1: Cold-induction time series...")

key_strains <- c("SW2B2","SW2C3","H2E-10C","SW1C-HB20","M2-6B2")

fig1_data <- master %>%
  filter(strain %in% key_strains,
         salinity_m==0,
         !is.na(delta_gamma_OD),
         !is.na(temperature)) %>%
  mutate(strain=factor(strain, levels=key_strains))

p_fig1 <- ggplot(fig1_data,
                 aes(x=day, y=delta_gamma_OD,
                     colour=temperature, group=temperature)) +
  geom_line(linewidth=0.8, alpha=0.9) +
  geom_point(size=1.8, alpha=0.8) +
  facet_wrap(~strain, nrow=2, scales="free_y") +
  scale_colour_manual(values=COL_TEMP, name="Temperature") +
  scale_x_continuous(breaks=c(1,3,5,7,9,11)) +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Day", y=expression(Delta*gamma*"/OD (mN/m per OD)"),
       title="Cold-induction of biosurfactant production",
       subtitle="Key cold-stimulated strains: temperature-only experiments") +
  THEME_PUB +
  theme(legend.position="bottom",
        strip.text=element_text(face="bold", size=9))

save_fig("Fig1_cold_induction_timeseries", p_fig1,
         w=FIG_W_DOUBLE, h=5.5, dir=OUT_MAIN)
message("  Fig 1 saved")

# =============================================================================
# FIGURE 2: Strain Classification — CSI ranking + producer classes
# =============================================================================
message("\nBuilding Figure 2: Strain classification...")

# Panel A: CSI lollipop chart
p_fig2a <- strain_class %>%
  filter(!is.na(CSI), !is.infinite(CSI)) %>%
  mutate(strain=fct_reorder(strain, CSI),
         CSI_class=factor(CSI_class,
                          levels=c("Cold-stimulated","Neutral","Warm-preferred"))) %>%
  ggplot(aes(x=CSI, y=strain, colour=CSI_class)) +
  geom_segment(aes(x=0, xend=CSI, y=strain, yend=strain),
               linewidth=0.6, alpha=0.7) +
  geom_point(size=3) +
  geom_vline(xintercept=c(0.7,1.5), linetype="dashed",
             colour="grey40", linewidth=0.5) +
  scale_colour_manual(values=COL_CSI, name="CSI class") +
  scale_x_continuous(labels=scales::comma) +
  labs(x="Cold Stimulation Index (AUC 4C / AUC 30C)",
       y=NULL,
       tag="A") +
  THEME_PUB +
  theme(legend.position="none")

# Panel B: Max production at 4C vs 30C scatter
p_fig2b <- strain_class %>%
  filter(!is.na(max_at_4C), !is.na(max_at_30C),
         !is.infinite(max_at_4C), !is.infinite(max_at_30C)) %>%
  ggplot(aes(x=max_at_30C, y=max_at_4C, colour=CSI_class, label=strain)) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              colour="grey50", linewidth=0.6) +
  geom_point(size=2.5, alpha=0.85) +
  geom_text_repel(size=2.2, max.overlaps=15,
                  segment.colour="grey70", segment.size=0.3) +
  scale_colour_manual(values=COL_CSI, name="CSI class", na.value="grey70") +
  scale_x_continuous(labels=scales::comma) +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Max(dg/OD) at 30C",
       y="Max(dg/OD) at 4C",
       tag="B") +
  THEME_PUB +
  theme(legend.position="right")

p_fig2 <- p_fig2a | p_fig2b
save_fig("Fig2_strain_classification", p_fig2,
         w=FIG_W_DOUBLE, h=5.5, dir=OUT_MAIN)
message("  Fig 2 saved")

# =============================================================================
# FIGURE 3: Salinity modulation — interaction plots
# =============================================================================
message("\nBuilding Figure 3: Salinity modulation...")

# Panel A: Interaction plot (emmeans) — refit LMM
df_sal <- summary_df %>%
  filter(exp_type=="temp_salinity",
         !is.na(max_delta_gamma_OD), !is.na(temperature)) %>%
  mutate(log_max_dg_OD=log10(max_delta_gamma_OD+1),
         salinity_m=factor(salinity_m, levels=c(0.5,1.0,1.5)),
         strain=factor(strain), experiment=factor(experiment))

lmm_sal <- lmer(log_max_dg_OD ~ temperature*salinity_m +
                  (1|strain)+(1|experiment),
                data=df_sal, REML=TRUE)

emm_grid <- emmeans(lmm_sal, ~temperature*salinity_m) %>%
  as.data.frame() %>%
  rename(temp=temperature, sal=salinity_m, estimate=emmean, se=SE) %>%
  mutate(sal=factor(as.character(sal), levels=c("0.5","1","1.5")),
         temp=factor(as.character(temp), levels=TEMP_LEVELS))

p_fig3a <- ggplot(emm_grid,
                  aes(x=temp, y=estimate,
                      colour=sal, group=sal)) +
  geom_line(linewidth=0.9, position=position_dodge(0.2)) +
  geom_point(size=3, position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se),
                width=0.12, linewidth=0.6,
                position=position_dodge(0.2)) +
  scale_colour_manual(values=COL_SAL, name="NaCl (M)") +
  labs(x="Temperature",
       y="log10(Max(dg/OD)+1)\nEMM ± SE",
       tag="A") +
  THEME_PUB

# Panel B: Bliss synergy plot
p_fig3b <- bliss %>%
  mutate(strain=fct_reorder(strain, synergy),
         syn_class=factor(syn_class,
                          levels=c("Synergistic","Additive","Antagonistic"))) %>%
  ggplot(aes(x=synergy, y=strain, fill=syn_class)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=c(-0.1,0.1), linetype="dashed",
             colour="grey40", linewidth=0.5) +
  scale_fill_manual(values=c(Synergistic="#1D9E75",
                              Additive="#888780",
                              Antagonistic="#D85A30"),
                    name="Interaction") +
  labs(x="Bliss synergy index\n(observed - expected)",
       y=NULL,
       tag="B") +
  THEME_PUB

# Panel C: Top producers at 4C + 1.5M
p_fig3c <- strain_prof %>%
  filter(temperature=="4C",
         as.character(salinity_m)=="1.5") %>%
  mutate(strain=fct_reorder(strain, mean_dg_OD)) %>%
  ggplot(aes(x=mean_dg_OD, y=strain, fill=CSI_class)) +
  geom_col(width=0.7) +
  scale_fill_manual(values=COL_CSI, name="CSI class",
                    na.value="grey70") +
  scale_x_continuous(labels=scales::comma) +
  labs(x="Mean Max(dg/OD)",
       y=NULL,
       tag="C",
       title="Top producers: 4C + 1.5M NaCl") +
  THEME_PUB

p_fig3 <- (p_fig3a | p_fig3b) / p_fig3c +
  plot_annotation(title="Salinity modulation of cold biosurfactant production")

save_fig("Fig3_salinity_modulation", p_fig3,
         w=FIG_W_DOUBLE, h=8, dir=OUT_MAIN)
message("  Fig 3 saved")

# =============================================================================
# FIGURE 4: Genomic integration — the genotype-phenotype bridge
# =============================================================================
message("\nBuilding Figure 4: Genomic integration...")

# Load WGS data for this figure
stress_genes <- read_tsv(FILE_STRESS_GENES, show_col_types=FALSE)
wgs_meta     <- read_tsv(FILE_WGS_META,     show_col_types=FALSE)

phys_4C <- summary_df %>%
  filter(!is.na(max_delta_gamma_OD), temperature=="4C") %>%
  group_by(strain) %>%
  summarise(max_dg_OD_4C=max(max_delta_gamma_OD,na.rm=TRUE),
            .groups="drop") %>%
  mutate(max_dg_OD_4C=if_else(is.infinite(max_dg_OD_4C),
                               NA_real_,max_dg_OD_4C))

gp_data <- phys_4C %>%
  left_join(meta %>% select(strain_phys,wgs_code,environment,
                             has_lipopeptide_bgc,bgc_type_lipopeptide,
                             has_wgs,CSI_class),
            by=c("strain"="strain_phys")) %>%
  filter(has_wgs==TRUE)

# Panel A: BGC presence vs production
p_fig4a <- gp_data %>%
  filter(!is.na(has_lipopeptide_bgc), !is.na(max_dg_OD_4C)) %>%
  mutate(bgc_label=if_else(has_lipopeptide_bgc,
                            "BGC+\n(Lichenysin/Surfactin)",
                            "BGC-")) %>%
  ggplot(aes(x=bgc_label, y=max_dg_OD_4C, colour=bgc_label)) +
  geom_boxplot(aes(fill=bgc_label), alpha=0.12, outlier.shape=NA,
               linewidth=0.5, width=0.45) +
  geom_beeswarm(size=2.5, cex=2.5, alpha=0.85) +
  geom_text_repel(aes(label=strain), size=2.0, max.overlaps=10,
                  segment.size=0.3, segment.colour="grey70") +
  scale_colour_manual(values=c("BGC+\n(Lichenysin/Surfactin)"="#1D9E75",
                                "BGC-"="#D85A30"), guide="none") +
  scale_fill_manual(values=c("BGC+\n(Lichenysin/Surfactin)"="#1D9E75",
                              "BGC-"="#D85A30"), guide="none") +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Lipopeptide BGC", y="Max(dg/OD) at 4C",
       tag="A") +
  THEME_PUB

# Panel B: Stress gene regulatory capacity
stress_summary <- stress_genes %>%
  group_by(Isolate,Category) %>%
  summarise(n=n(),.groups="drop") %>%
  pivot_wider(names_from=Category,values_from=n,values_fill=0)

bs_col <- names(stress_summary)[str_detect(
  names(stress_summary), regex("biosurfactant",ignore_case=TRUE))]

if (length(bs_col)==1) {
  stress_phys <- stress_summary %>%
    left_join(wgs_meta %>% select(WGS_Code,Strain_Code),
              by=c("Isolate"="WGS_Code")) %>%
    inner_join(phys_4C, by=c("Strain_Code"="strain")) %>%
    filter(!is.na(max_dg_OD_4C)) %>%
    left_join(meta %>% select(strain_phys,CSI_class),
              by=c("Strain_Code"="strain_phys"))

  p_fig4b <- ggplot(stress_phys,
                    aes(x=.data[[bs_col]], y=max_dg_OD_4C,
                        colour=CSI_class, label=Strain_Code)) +
    geom_smooth(method="lm", se=TRUE, colour="grey50",
                linewidth=0.6, alpha=0.15) +
    geom_point(size=2.5, alpha=0.85) +
    geom_text_repel(size=2.0, max.overlaps=10,
                    segment.size=0.3, segment.colour="grey70") +
    scale_colour_manual(values=COL_CSI, name="CSI class",
                        na.value="grey70") +
    scale_y_continuous(labels=scales::comma) +
    annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5,
             label=paste0("rho=0.874\np=0.0002"),
             size=3, colour="grey30") +
    labs(x="Biosurfactant regulation\ngene count",
         y="Max(dg/OD) at 4C",
         tag="B") +
    THEME_PUB
} else {
  p_fig4b <- ggplot() + labs(title="Stress gene data unavailable") + THEME_PUB
}

# Panel C: Environment comparison
p_fig4c <- gp_data %>%
  filter(!is.na(environment), !is.na(max_dg_OD_4C)) %>%
  ggplot(aes(x=environment, y=max_dg_OD_4C, colour=environment)) +
  geom_boxplot(aes(fill=environment), alpha=0.12, outlier.shape=NA,
               linewidth=0.5, width=0.45) +
  geom_beeswarm(size=2.5, cex=2.5, alpha=0.85) +
  geom_text_repel(aes(label=strain), size=2.0, max.overlaps=8,
                  segment.size=0.3, segment.colour="grey70") +
  scale_colour_manual(values=COL_ENV, guide="none") +
  scale_fill_manual(values=COL_ENV,   guide="none") +
  scale_y_continuous(labels=scales::comma) +
  annotate("text", x=1.5, y=max(gp_data$max_dg_OD_4C,na.rm=TRUE)*0.9,
           label="p=0.038", size=3, colour="grey30") +
  labs(x="Ecological origin", y="Max(dg/OD) at 4C",
       tag="C") +
  THEME_PUB

# Panel D: CSI class by BGC status (stacked bar)
p_fig4d <- gp_data %>%
  filter(!is.na(has_lipopeptide_bgc), !is.na(CSI_class)) %>%
  mutate(bgc_label=if_else(has_lipopeptide_bgc,"BGC+","BGC-"),
         CSI_class=factor(CSI_class,
                          levels=c("Cold-stimulated","Neutral","Warm-preferred"))) %>%
  count(bgc_label, CSI_class) %>%
  ggplot(aes(x=bgc_label, y=n, fill=CSI_class)) +
  geom_col(position="fill", width=0.6) +
  scale_fill_manual(values=COL_CSI, name="CSI class") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Lipopeptide BGC", y="Proportion",
       tag="D") +
  THEME_PUB

p_fig4 <- (p_fig4a | p_fig4b) / (p_fig4c | p_fig4d) +
  plot_annotation(
    title="Genotype-phenotype integration",
    subtitle="Linking BGC presence and regulatory capacity to cold biosurfactant production"
  )

save_fig("Fig4_genomic_integration", p_fig4,
         w=FIG_W_DOUBLE, h=8, dir=OUT_MAIN)
message("  Fig 4 saved")

# =============================================================================
# FIGURE 5: Ecological filter — full strain heatmap
# =============================================================================
message("\nBuilding Figure 5: Ecological filter heatmap...")

hmat_data <- summary_df %>%
  filter(!is.na(temperature), !is.na(max_delta_gamma_OD),
         salinity_m==0) %>%
  group_by(strain, temperature) %>%
  summarise(val=mean(max_delta_gamma_OD,na.rm=TRUE),.groups="drop") %>%
  pivot_wider(names_from=temperature, values_from=val, values_fill=0) %>%
  column_to_rownames("strain")

col_order <- intersect(TEMP_LEVELS, names(hmat_data))
hmat_data <- hmat_data[,col_order,drop=FALSE]

ann_row <- strain_class %>%
  select(strain, producer_class, CSI_class, environment) %>%
  filter(strain %in% rownames(hmat_data)) %>%
  distinct(strain,.keep_all=TRUE) %>%
  column_to_rownames("strain")

ann_colors <- list(
  producer_class = c("High"="#1D9E75","Moderate"="#EF9F27",
                     "Low"="#D3D1C7","Non-producer"="grey90"),
  CSI_class      = COL_CSI,
  environment    = COL_ENV
)

hmat_matrix <- as.matrix(hmat_data)
hmat_matrix[is.na(hmat_matrix)] <- 0

out_fig5 <- file.path(OUT_MAIN,"Fig5_ecological_filter_heatmap.pdf")
pdf(out_fig5, width=FIG_W_SINGLE*1.4, height=FIG_H_MAX*0.75)
pheatmap(
  hmat_matrix,
  annotation_row    = ann_row[rownames(hmat_matrix),,drop=FALSE],
  annotation_colors = ann_colors,
  color             = colorRampPalette(
    c("white","#E8F5F0","#9FE1CB","#1D9E75","#085041"))(60),
  clustering_method = "ward.D2",
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize          = 9,
  fontsize_row      = 8,
  angle_col         = 0,
  main              = "Biosurfactant production across temperatures — all 22 strains",
  na_col            = "grey95",
  border_color      = "white"
)
dev.off()
message("  Fig 5 saved: ", out_fig5)

# =============================================================================
# SUPPLEMENTARY FIGURES
# =============================================================================
message("\nBuilding supplementary figures...")

# S1: Baseline ST homogeneity
baseline_data <- master %>%
  filter(day==1) %>%
  mutate(experiment=factor(experiment,
         levels=paste0("exp",2:11)))

p_s1 <- ggplot(baseline_data, aes(x=experiment, y=st)) +
  geom_boxplot(fill="grey90",colour="grey40",outlier.shape=NA,
               linewidth=0.4,width=0.5) +
  geom_jitter(aes(colour=temperature),width=0.15,size=1.5,alpha=0.7) +
  geom_hline(yintercept=c(68,71),linetype="dotted",
             colour="red",linewidth=0.5) +
  scale_colour_manual(values=COL_TEMP,name="Temperature") +
  scale_y_continuous(limits=c(67,72),breaks=seq(67,72,1)) +
  labs(x="Experiment",y="Baseline ST (mN/m)",
       title="S1: Baseline surface tension homogeneity across experiments") +
  THEME_PUB +
  theme(axis.text.x=element_text(angle=45,hjust=1))

save_fig("FigS1_baseline_homogeneity", p_s1,
         w=FIG_W_DOUBLE, h=3.5, dir=OUT_SUPP)

# S2: GMM producer classification
gmm_thresholds <- readRDS(
  file.path(RESULTS,"module02_exploratory","gmm_thresholds.rds"))

all_vals <- summary_df %>%
  filter(!is.na(max_delta_gamma_OD), max_delta_gamma_OD>0) %>%
  pull(max_delta_gamma_OD)

p_s2 <- ggplot(data.frame(x=all_vals), aes(x=x)) +
  geom_histogram(aes(y=after_stat(density)), bins=40,
                 fill="grey80",colour="grey60",linewidth=0.3) +
  geom_density(colour="grey40",linewidth=0.6) +
  geom_vline(xintercept=gmm_thresholds$thresh_low_moderate,
             linetype="dashed",colour="#EF9F27",linewidth=0.8) +
  geom_vline(xintercept=gmm_thresholds$thresh_moderate_high,
             linetype="dashed",colour="#1D9E75",linewidth=0.8) +
  annotate("text",x=gmm_thresholds$thresh_low_moderate*0.4,
           y=Inf,vjust=1.5,size=3,label="Low",colour="#EF9F27") +
  annotate("text",x=(gmm_thresholds$thresh_low_moderate+
                       gmm_thresholds$thresh_moderate_high)/2,
           y=Inf,vjust=1.5,size=3,label="Moderate",colour="grey40") +
  annotate("text",x=gmm_thresholds$thresh_moderate_high*1.4,
           y=Inf,vjust=1.5,size=3,label="High",colour="#1D9E75") +
  labs(x=expression("Max("*Delta*gamma*"/OD)"),
       y="Density",
       title="S2: Gaussian mixture model — producer classification") +
  THEME_PUB

save_fig("FigS2_GMM_classification", p_s2,
         w=FIG_W_SINGLE, h=3.5, dir=OUT_SUPP)

# S3: Full time-series all strains (temp-only)
p_s3 <- master %>%
  filter(salinity_m==0, !is.na(delta_gamma_OD), !is.na(temperature)) %>%
  ggplot(aes(x=day, y=delta_gamma_OD,
             colour=temperature, group=interaction(temperature,experiment))) +
  geom_line(linewidth=0.5, alpha=0.7) +
  facet_wrap(~strain, ncol=5, scales="free_y") +
  scale_colour_manual(values=COL_TEMP, name="Temperature") +
  scale_x_continuous(breaks=c(1,5,9,11)) +
  labs(x="Day", y=expression(Delta*gamma*"/OD"),
       title="S3: Full time-series — all strains, temperature-only experiments") +
  THEME_PUB +
  theme(legend.position="bottom",
        strip.text=element_text(size=6),
        axis.text=element_text(size=7))

save_fig("FigS3_full_timeseries", p_s3,
         w=FIG_W_DOUBLE, h=FIG_H_MAX*0.8, dir=OUT_SUPP)

# S4: OD vs production decoupling
p_s4 <- summary_df %>%
  filter(salinity_m==0, !is.na(max_OD), !is.na(max_delta_gamma_OD),
         !is.na(temperature)) %>%
  ggplot(aes(x=max_OD, y=max_delta_gamma_OD,
             colour=temperature, label=strain)) +
  geom_smooth(method="lm", se=TRUE, linewidth=0.6,
              alpha=0.12, colour="grey50") +
  geom_point(size=2, alpha=0.85) +
  geom_text_repel(size=2.0, max.overlaps=8,
                  segment.size=0.3, segment.colour="grey70") +
  scale_colour_manual(values=COL_TEMP, name="Temperature") +
  scale_y_continuous(labels=scales::comma) +
  facet_wrap(~temperature, nrow=1) +
  labs(x="Max OD", y=expression("Max("*Delta*gamma*"/OD)"),
       title="S4: Growth vs biosurfactant production — decoupling analysis") +
  THEME_PUB + theme(legend.position="none")

save_fig("FigS4_OD_vs_production", p_s4,
         w=FIG_W_DOUBLE, h=3.8, dir=OUT_SUPP)

# =============================================================================
# Summary
# =============================================================================
message("\n=== MODULE 7 SUMMARY ===")
message("Main figures (", OUT_MAIN, "):")
main_figs <- list.files(OUT_MAIN, pattern="Fig[0-9].*\\.pdf")
for(f in main_figs) message("  ",f)
message("\nSupplementary figures (", OUT_SUPP, "):")
supp_figs <- list.files(OUT_SUPP, pattern=".*\\.pdf")
for(f in supp_figs) message("  ",f)
message("\n=== MODULE 7 COMPLETE ===\n")
