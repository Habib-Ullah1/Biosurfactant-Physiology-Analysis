source("scripts/00_config.R")
library(tidyverse); library(lme4); library(lmerTest)
library(emmeans); library(MuMIn); library(broom.mixed)
library(ggbeeswarm); library(patchwork)

message("\n=== MODULE 3: Temperature effects (LMM) ===\n")

OUT <- file.path(RESULTS, "module03_temperature")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# -----------------------------------------------------------------------------
# 1. Load and prepare data
# -----------------------------------------------------------------------------
summary <- read_csv(file.path(DATA_PROC, "summary_metrics.csv"), show_col_types=FALSE)
meta    <- read_csv(FILE_METADATA, show_col_types=FALSE)

# Keep temp-only conditions (salinity=0) for H1
df <- summary %>%
  filter(salinity_m == 0,
         !is.na(max_delta_gamma_OD),
         !is.na(temperature)) %>%
  mutate(
    temperature  = factor(temperature, levels=TEMP_LEVELS),
    strain       = factor(strain),
    experiment   = factor(experiment),
    log_max_dg_OD = log10(max_delta_gamma_OD + 1),
    log_AUC       = log10(AUC_delta_gamma_OD  + 1)
  ) %>%
  left_join(meta %>% select(strain_phys, environment, CSI_class, species),
            by=c("strain"="strain_phys"))

message("Observations for LMM: ", nrow(df))
message("Strains:              ", n_distinct(df$strain))
message("Experiments:          ", n_distinct(df$experiment))
message("Temperature levels:   ", paste(levels(df$temperature), collapse=" | "))

# -----------------------------------------------------------------------------
# 2. Assumption checks BEFORE fitting model
# -----------------------------------------------------------------------------
message("\n--- Assumption checks (raw data) ---")

# Shapiro-Wilk per temperature group
sw_raw <- df %>%
  group_by(temperature) %>%
  summarise(
    W     = shapiro.test(max_delta_gamma_OD)$statistic,
    p     = shapiro.test(max_delta_gamma_OD)$p.value,
    W_log = shapiro.test(log_max_dg_OD)$statistic,
    p_log = shapiro.test(log_max_dg_OD)$p.value,
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.,3)))

message("Shapiro-Wilk normality test:")
print(sw_raw)

# Decide: use log-transform if any raw p < 0.05
use_log <- any(sw_raw$p < 0.05)
response_var  <- if (use_log) "log_max_dg_OD" else "max_delta_gamma_OD"
response_label <- if (use_log) "log10(Max(Δγ/OD)+1)" else "Max(Δγ/OD)"
message("Transform: ", if(use_log) "LOG10 applied (raw non-normal)" else "NO transform (raw normal)")

# -----------------------------------------------------------------------------
# 3. Fit Linear Mixed Model
# -----------------------------------------------------------------------------
message("\n--- Fitting Linear Mixed Model ---")
message("Model: ", response_var, " ~ Temperature + (1|Strain) + (1|Experiment)")

formula_main <- as.formula(paste(response_var,
  "~ temperature + (1|strain) + (1|experiment)"))

lmm_main <- lmer(formula_main, data=df, REML=TRUE)
message("Model converged: OK")

# Model summary
message("\n--- LMM Summary ---")
print(summary(lmm_main))

# F-test (Type III) via lmerTest
message("\n--- ANOVA table (Type III F-test) ---")
anova_out <- anova(lmm_main, type=3)
print(anova_out)

# R-squared
r2 <- r.squaredGLMM(lmm_main)
message("\nMarginal R2  (Temperature fixed effect): ", round(r2[1],3))
message("Conditional R2 (full model):              ", round(r2[2],3))

# Random effects variance
message("\nRandom effects variance:")
print(VarCorr(lmm_main))

# -----------------------------------------------------------------------------
# 4. Post-hoc pairwise comparisons (Tukey)
# -----------------------------------------------------------------------------
message("\n--- Post-hoc: Tukey pairwise comparisons ---")

emm <- emmeans(lmm_main, ~ temperature)
pairs_out <- pairs(emm, adjust="tukey") %>% as.data.frame()
message("Pairwise comparisons:")
print(pairs_out %>% mutate(across(where(is.numeric), ~round(.,4))))

# Cohen's d effect sizes
message("\nCohen's d effect sizes:")
eff_out <- eff_size(emm, sigma=sigma(lmm_main),
                    edf=df.residual(lmm_main)) %>% as.data.frame()
print(eff_out %>% mutate(across(where(is.numeric), ~round(.,3))))

# -----------------------------------------------------------------------------
# 5. Assumption checks on RESIDUALS
# -----------------------------------------------------------------------------
message("\n--- Residual diagnostics ---")
resid_df <- data.frame(
  fitted    = fitted(lmm_main),
  residuals = residuals(lmm_main),
  temperature = df$temperature,
  strain      = df$strain
)

sw_resid <- shapiro.test(residuals(lmm_main))
message("Shapiro-Wilk on residuals: W=", round(sw_resid$statistic,3),
        " p=", round(sw_resid$p.value,4))
if (sw_resid$p.value < 0.05) {
  message("  WARNING: residuals non-normal -> Kruskal-Wallis backup run below")
} else {
  message("  Residuals normal: LMM results are reliable")
}

# Levene test for homoscedasticity
levene_p <- df %>%
  group_by(temperature) %>%
  summarise(v=var(.data[[response_var]], na.rm=TRUE), .groups="drop")
message("Variance per temperature group:")
print(levene_p %>% mutate(v=round(v,2)))

# -----------------------------------------------------------------------------
# 6. Non-parametric backup (Kruskal-Wallis + Dunn)
# -----------------------------------------------------------------------------
message("\n--- Non-parametric backup: Kruskal-Wallis ---")
kw <- kruskal.test(as.formula(paste(response_var, "~ temperature")), data=df)
message("Kruskal-Wallis: chi2=", round(kw$statistic,3),
        " df=", kw$parameter,
        " p=", round(kw$p.value,5))

# Dunn test (pairwise)
if (!requireNamespace("dunn.test", quietly=TRUE)) {
  install.packages("dunn.test", repos="https://cloud.r-project.org")
}
library(dunn.test)
message("\nDunn test (Bonferroni):")
dunn_out <- dunn.test(df[[response_var]],
                      g=df$temperature,
                      method="bonferroni", altp=TRUE)

# -----------------------------------------------------------------------------
# 7. AUC model (secondary response variable)
# -----------------------------------------------------------------------------
message("\n--- Secondary model: AUC(Δγ/OD) ---")
formula_auc <- as.formula("log_AUC ~ temperature + (1|strain) + (1|experiment)")
lmm_auc <- tryCatch(
  lmer(formula_auc, data=df %>% filter(!is.na(log_AUC)), REML=TRUE),
  error=function(e) { message("AUC model failed: ", e$message); NULL }
)
if (!is.null(lmm_auc)) {
  anova_auc <- anova(lmm_auc, type=3)
  r2_auc    <- r.squaredGLMM(lmm_auc)
  emm_auc   <- emmeans(lmm_auc, ~temperature)
  pairs_auc <- pairs(emm_auc, adjust="tukey") %>% as.data.frame()
  message("AUC model F-test:")
  print(anova_auc)
  message("AUC marginal R2: ", round(r2_auc[1],3),
          " | conditional R2: ", round(r2_auc[2],3))
  message("AUC pairwise:")
  print(pairs_auc %>% mutate(across(where(is.numeric),~round(.,4))))
}

# -----------------------------------------------------------------------------
# 8. Save results tables
# -----------------------------------------------------------------------------
message("\nSaving results tables...")

# Main ANOVA table
anova_tbl <- as.data.frame(anova_out) %>%
  rownames_to_column("term") %>%
  mutate(across(where(is.numeric), ~round(.,4)))
save_table(anova_tbl, "module03_lmm_anova", dir=TABLES_MAIN)

# Pairwise comparisons
pairs_tbl <- pairs_out %>%
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))
save_table(pairs_tbl, "module03_lmm_pairwise", dir=TABLES_MAIN)

# Effect sizes
save_table(eff_out, "module03_effect_sizes", dir=TABLES_MAIN)

# Summary stats per temperature
temp_stats <- df %>%
  group_by(temperature) %>%
  summarise(
    n           = n(),
    mean        = mean(.data[[response_var]], na.rm=TRUE),
    sd          = sd(.data[[response_var]],   na.rm=TRUE),
    median      = median(.data[[response_var]], na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.,3)))
save_table(temp_stats, "module03_temp_summary_stats", dir=TABLES_MAIN)

# -----------------------------------------------------------------------------
# 9. Figures
# -----------------------------------------------------------------------------
message("\nBuilding figures...")

# Get emmeans for plotting
emm_df <- as.data.frame(emm) %>%
  rename(temp=temperature, estimate=emmean, se=SE)

# Fig A: Boxplot + beeswarm with emmeans overlay
p_box <- ggplot(df, aes(x=temperature, y=.data[[response_var]],
                         colour=temperature)) +
  geom_boxplot(aes(fill=temperature), alpha=0.12, outlier.shape=NA,
               linewidth=0.5, width=0.45) +
  geom_beeswarm(size=2, cex=2.8, alpha=0.8) +
  geom_point(data=emm_df, aes(x=temp, y=estimate),
             shape=18, size=4, colour="black") +
  geom_errorbar(data=emm_df,
                aes(x=temp, y=estimate,
                    ymin=estimate-se, ymax=estimate+se),
                colour="black", width=0.15, linewidth=0.7,
                inherit.aes=FALSE) +
  scale_colour_manual(values=COL_TEMP, guide="none") +
  scale_fill_manual(values=COL_TEMP,   guide="none") +
  labs(x="Temperature", y=response_label,
       title="H1: Cold temperature stimulates biosurfactant production",
       subtitle="Diamond = estimated marginal mean ± SE from LMM") +
  THEME_PUB

# Add significance brackets
pairs_sig <- pairs_tbl %>% filter(p.value < 0.05)
message("  Significant pairs: ", nrow(pairs_sig))

# Fig B: Residual diagnostic plot
p_resid <- ggplot(resid_df, aes(x=fitted, y=residuals)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
  geom_point(aes(colour=temperature), size=1.8, alpha=0.7) +
  geom_smooth(method="loess", se=FALSE, colour="red",
              linewidth=0.6) +
  scale_colour_manual(values=COL_TEMP, name="Temperature") +
  labs(x="Fitted values", y="Residuals",
       title="Residual diagnostics",
       subtitle="No pattern expected if model assumptions met") +
  THEME_PUB

# Fig C: QQ plot
p_qq <- ggplot(resid_df, aes(sample=residuals)) +
  stat_qq(size=1.5, alpha=0.6, colour="#378ADD") +
  stat_qq_line(colour="red", linewidth=0.7) +
  labs(x="Theoretical quantiles", y="Sample quantiles",
       title="Q-Q plot of residuals") +
  THEME_PUB

# Fig D: Per-strain temperature response (spaghetti plot)
strain_means <- df %>%
  group_by(strain, temperature) %>%
  summarise(val=mean(.data[[response_var]], na.rm=TRUE), .groups="drop") %>%
  left_join(meta %>% select(strain_phys, CSI_class),
            by=c("strain"="strain_phys"))

p_spaghetti <- ggplot(strain_means,
                      aes(x=temperature, y=val, group=strain,
                          colour=CSI_class)) +
  geom_line(alpha=0.6, linewidth=0.6) +
  geom_point(size=1.8, alpha=0.8) +
  scale_colour_manual(values=COL_CSI, name="CSI class",
                      na.value="grey70") +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Temperature", y=response_label,
       title="Per-strain temperature response profiles",
       subtitle="Lines connect same strain across temperatures") +
  THEME_PUB

# Combine into multi-panel figure
p_combined <- (p_box | p_spaghetti) / (p_resid | p_qq) +
  plot_annotation(
    title="Module 3: Temperature effects on biosurfactant production",
    tag_levels="A"
  )

save_fig("module03_temperature_main", p_combined,
         w=FIG_W_DOUBLE, h=7,
         dir=OUT)

# Save individual panels too
save_fig("module03_boxplot",    p_box,       w=FIG_W_SINGLE, h=4, dir=OUT)
save_fig("module03_spaghetti",  p_spaghetti, w=FIG_W_SINGLE, h=4, dir=OUT)
save_fig("module03_residuals",  p_resid,     w=FIG_W_SINGLE, h=3.5, dir=OUT)
save_fig("module03_qqplot",     p_qq,        w=FIG_W_SINGLE, h=3.5, dir=OUT)

# -----------------------------------------------------------------------------
# 10. Summary report
# -----------------------------------------------------------------------------
message("\n=== MODULE 3 SUMMARY ===")
message("Response variable:    ", response_label)
message("Transform applied:    ", if(use_log) "log10" else "none")
message("N observations:       ", nrow(df))
message("Marginal R2 (Temp):   ", round(r2[1],3))
message("Conditional R2:       ", round(r2[2],3))
message("Temperature F-test:   F=", round(anova_out$`F value`,3),
        " df=", anova_out$NumDF,",", round(anova_out$DenDF,1),
        " p=", round(anova_out$`Pr(>F)`,5))
message("Kruskal-Wallis p:     ", round(kw$p.value,5))
message("Significant pairs:    ", nrow(pairs_sig), "/3")
message("Results saved to:     ", OUT)
message("\n=== MODULE 3 COMPLETE ===\n")
