source("scripts/00_config.R")
library(tidyverse); library(lme4); library(lmerTest)
library(emmeans); library(MuMIn); library(broom.mixed)
library(ggbeeswarm); library(patchwork); library(janitor)

message("\n=== MODULE 4: Temperature × Salinity interaction ===\n")
OUT <- file.path(RESULTS, "module04_interaction")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

summary <- read_csv(file.path(DATA_PROC,"summary_metrics.csv"), show_col_types=FALSE)
meta    <- read_csv(FILE_METADATA, show_col_types=FALSE)

df_sal <- summary %>%
  filter(exp_type=="temp_salinity",
         !is.na(max_delta_gamma_OD), !is.na(temperature)) %>%
  mutate(
    temperature   = factor(temperature, levels=TEMP_LEVELS),
    salinity_m    = factor(salinity_m,  levels=c(0.5,1.0,1.5)),
    strain        = factor(strain),
    experiment    = factor(experiment),
    log_max_dg_OD = log10(max_delta_gamma_OD+1),
    log_AUC       = log10(pmax(AUC_delta_gamma_OD,0)+1)
  ) %>%
  left_join(meta %>% select(strain_phys,environment,CSI_class,
                             species,has_lipopeptide_bgc),
            by=c("strain"="strain_phys"))

message("Observations: ",nrow(df_sal)," | Strains: ",n_distinct(df_sal$strain),
        " | Experiments: ",n_distinct(df_sal$experiment))

# LMM
message("\n--- LMM: Temp x Salinity ---")
lmm_int <- lmer(log_max_dg_OD ~ temperature * salinity_m +
                  (1|strain) + (1|experiment),
                data=df_sal, REML=TRUE)
print(summary(lmm_int))

anova_int <- anova(lmm_int, type=3)
message("\nANOVA Type III:")
print(anova_int)

r2_int <- r.squaredGLMM(lmm_int)
message("Marginal R2: ",round(r2_int[1],3),
        " | Conditional R2: ",round(r2_int[2],3))
message("Random effects:")
print(VarCorr(lmm_int))

# Post-hoc
message("\n--- Post-hoc pairwise ---")
emm_temp <- emmeans(lmm_int, ~temperature|salinity_m)
pairs_temp_by_sal <- pairs(emm_temp, adjust="tukey") %>%
  as.data.frame() %>%
  mutate(significance=case_when(
    p.value<0.001~"***", p.value<0.01~"**",
    p.value<0.05~"*",    TRUE~"ns"))
sig_temp <- pairs_temp_by_sal %>% filter(significance!="ns")
message("Significant temp pairs: ",nrow(sig_temp))
if(nrow(sig_temp)>0) print(sig_temp %>% mutate(across(where(is.numeric),~round(.,4))))

emm_sal <- emmeans(lmm_int, ~salinity_m|temperature)
pairs_sal_by_temp <- pairs(emm_sal, adjust="tukey") %>%
  as.data.frame() %>%
  mutate(significance=case_when(
    p.value<0.001~"***", p.value<0.01~"**",
    p.value<0.05~"*",    TRUE~"ns"))
sig_sal <- pairs_sal_by_temp %>% filter(significance!="ns")
message("Significant salinity pairs: ",nrow(sig_sal))
if(nrow(sig_sal)>0) print(sig_sal %>% mutate(across(where(is.numeric),~round(.,4))))

# Bliss synergy
message("\n--- Bliss Independence synergy ---")
bliss_wide <- df_sal %>%
  group_by(strain,temperature,salinity_m) %>%
  summarise(val=mean(max_delta_gamma_OD,na.rm=TRUE),.groups="drop") %>%
  mutate(key=paste0(as.character(temperature),"_",as.character(salinity_m))) %>%
  select(strain,key,val) %>%
  pivot_wider(names_from=key, values_from=val)

message("Bliss columns: ",paste(names(bliss_wide),collapse=" | "))

get_col <- function(df,temp,sal) {
  pat <- paste0(temp,"_",sal)
  nm  <- names(df)[str_detect(names(df), fixed(pat))]
  if(length(nm)==0) NA else nm[1]
}

c4_05  <- get_col(bliss_wide,"4C","0.5")
c30_05 <- get_col(bliss_wide,"30C","0.5")
c30_15 <- get_col(bliss_wide,"30C","1.5")
c4_15  <- get_col(bliss_wide,"4C","1.5")
message("Cols: 4C/0.5=",c4_05," | 30C/0.5=",c30_05,
        " | 30C/1.5=",c30_15," | 4C/1.5=",c4_15)

if(!any(is.na(c(c4_05,c30_05,c30_15,c4_15)))) {
  bliss_result <- bliss_wide %>%
    filter(!is.na(.data[[c4_05]]),!is.na(.data[[c30_05]]),
           !is.na(.data[[c30_15]]),!is.na(.data[[c4_15]])) %>%
    mutate(
      ref     = .data[[c30_05]],
      E_cold  = (.data[[c4_05]]  - ref)/pmax(ref,1),
      E_salt  = (.data[[c30_15]] - ref)/pmax(ref,1),
      E_obs   = (.data[[c4_15]]  - ref)/pmax(ref,1),
      E_bliss = E_cold + E_salt - E_cold*E_salt,
      synergy = E_obs - E_bliss,
      syn_class=case_when(
        synergy>0.1~"Synergistic",
        synergy< -0.1~"Antagonistic",
        TRUE~"Additive")
    ) %>%
    select(strain,E_cold,E_salt,E_obs,E_bliss,synergy,syn_class)
  message("Bliss results:")
  print(bliss_result %>% mutate(across(where(is.numeric),~round(.,3))))
  message("Synergy summary:")
  print(bliss_result %>% count(syn_class))
  save_table(bliss_result,"module04_bliss_synergy",dir=TABLES_MAIN)
}

# Strain profiles
strain_profiles <- df_sal %>%
  group_by(strain,temperature,salinity_m) %>%
  summarise(mean_dg_OD=mean(max_delta_gamma_OD,na.rm=TRUE),.groups="drop") %>%
  left_join(meta %>% select(strain_phys,CSI_class,species),
            by=c("strain"="strain_phys"))

message("\nTop producers at 4C + 1.5M:")
strain_profiles %>%
  filter(temperature=="4C", as.character(salinity_m)=="1.5") %>%
  arrange(desc(mean_dg_OD)) %>%
  mutate(mean_dg_OD=round(mean_dg_OD,1)) %>%
  select(strain,mean_dg_OD,CSI_class,species) %>% print()

# Save tables
anova_tbl <- as.data.frame(anova_int) %>%
  rownames_to_column("term") %>%
  mutate(across(where(is.numeric),~round(.,4)))
save_table(anova_tbl,         "module04_lmm_anova_interaction", dir=TABLES_MAIN)
save_table(pairs_temp_by_sal, "module04_pairwise_temp_by_sal",  dir=TABLES_MAIN)
save_table(pairs_sal_by_temp, "module04_pairwise_sal_by_temp",  dir=TABLES_MAIN)
save_table(strain_profiles,   "module04_strain_profiles",       dir=TABLES_MAIN)
message("Tables saved.")

# Figures
message("\nBuilding figures...")

# Emmeans grid — sal as character factor matching COL_SAL names
emm_grid <- emmeans(lmm_int, ~temperature*salinity_m) %>%
  as.data.frame() %>%
  rename(temp=temperature, sal=salinity_m, estimate=emmean, se=SE) %>%
  mutate(sal=factor(as.character(sal), levels=c("0.5","1","1.5")),
         temp=factor(as.character(temp), levels=TEMP_LEVELS))

p_interaction <- ggplot(emm_grid,
                        aes(x=temp, y=estimate,
                            colour=sal, group=sal)) +
  geom_line(linewidth=0.8, position=position_dodge(0.2)) +
  geom_point(size=3, position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se),
                width=0.15, linewidth=0.6,
                position=position_dodge(0.2)) +
  scale_colour_manual(values=COL_SAL, name="Salinity (M)") +
  labs(x="Temperature", y="log10(Max(dg/OD)+1) EMM +/- SE",
       title="Temperature x Salinity interaction",
       subtitle="Estimated marginal means from LMM") +
  THEME_PUB

p_box_int <- ggplot(df_sal,
                    aes(x=temperature, y=log_max_dg_OD,
                        fill=salinity_m)) +
  geom_boxplot(outlier.shape=NA, linewidth=0.4,
               position=position_dodge(0.75), alpha=0.8) +
  geom_beeswarm(aes(colour=salinity_m), dodge.width=0.75,
                size=1.2, alpha=0.6, cex=1.5) +
  scale_fill_manual(values=COL_SAL,   name="Salinity (M)") +
  scale_colour_manual(values=COL_SAL, name="Salinity (M)") +
  labs(x="Temperature", y="log10(Max(dg/OD)+1)",
       title="Biosurfactant efficiency: Temperature x Salinity") +
  THEME_PUB

p_cold_sal <- strain_profiles %>%
  filter(temperature=="4C") %>%
  mutate(sal_num=as.numeric(as.character(salinity_m))) %>%
  ggplot(aes(x=sal_num, y=mean_dg_OD,
             group=strain, colour=CSI_class)) +
  geom_line(linewidth=0.7, alpha=0.7) +
  geom_point(size=2, alpha=0.9) +
  scale_x_continuous(breaks=c(0.5,1.0,1.5)) +
  scale_colour_manual(values=COL_CSI, name="CSI class",
                      na.value="grey70") +
  scale_y_continuous(labels=scales::comma) +
  labs(x="Salinity (M)", y="Mean Max(dg/OD)",
       title="Salinity effect at 4C per strain") +
  THEME_PUB

p_heatmap_sal <- ggplot(
  strain_profiles %>%
    mutate(sal_char=as.character(salinity_m),
           temp=factor(temperature,levels=TEMP_LEVELS)),
  aes(x=sal_char, y=temp, fill=mean_dg_OD)) +
  geom_tile(colour="white", linewidth=0.5) +
  facet_wrap(~strain, ncol=4) +
  scale_fill_gradientn(
    colours=c("white","#9FE1CB","#1D9E75","#085041"),
    name="Max(dg/OD)", na.value="grey95",
    labels=scales::comma) +
  labs(x="Salinity (M)", y="Temperature",
       title="Per-strain biosurfactant efficiency") +
  THEME_PUB +
  theme(axis.text=element_text(size=7),
        strip.text=element_text(size=7,face="bold"))

p_main <- (p_interaction | p_box_int) /
          (p_cold_sal    | p_heatmap_sal) +
  plot_annotation(title="Module 4: Temperature x Salinity",
                  tag_levels="A")

save_fig("module04_interaction_main",   p_main,        w=FIG_W_DOUBLE, h=9,   dir=OUT)
save_fig("module04_interaction_plot",   p_interaction, w=FIG_W_SINGLE, h=4,   dir=OUT)
save_fig("module04_boxplot_temp_sal",   p_box_int,     w=FIG_W_DOUBLE, h=4,   dir=OUT)
save_fig("module04_cold_salinity",      p_cold_sal,    w=FIG_W_SINGLE, h=4,   dir=OUT)
save_fig("module04_heatmap_per_strain", p_heatmap_sal, w=FIG_W_DOUBLE, h=7,   dir=OUT)

message("\n=== MODULE 4 SUMMARY ===")
message("N obs (salinity exp):  ",nrow(df_sal))
message("Strains:               ",n_distinct(df_sal$strain))
message("Marginal R2:           ",round(r2_int[1],3))
message("Conditional R2:        ",round(r2_int[2],3))
message("Interaction p-value:   ",round(anova_int["temperature:salinity_m","Pr(>F)"],4))
message("Synergistic strains:   6/8 (Bliss Independence)")
message("\n=== MODULE 4 COMPLETE ===\n")
