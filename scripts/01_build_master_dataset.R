source("scripts/00_config.R")
library(tidyverse); library(readxl); library(janitor)
message("\n=== MODULE 1: Building master dataset ===\n")

day_map <- c(`1`=1,`2`=3,`3`=5,`4`=7,`5`=9,`6`=11)

experiments <- list(
  list(id="exp2",  file=FILE_TEMP, sheet="exp2 raw data",  type="temp_only",     has_salinity=FALSE),
  list(id="exp3",  file=FILE_TEMP, sheet="exp3 raw data",  type="temp_only",     has_salinity=FALSE),
  list(id="exp4",  file=FILE_TEMP, sheet="exp4 raw data",  type="temp_only",     has_salinity=FALSE),
  list(id="exp5",  file=FILE_TEMP, sheet="exp5 raw data",  type="temp_only",     has_salinity=FALSE),
  list(id="exp6",  file=FILE_TEMP, sheet="exp6 raw data",  type="temp_only",     has_salinity=FALSE),
  list(id="exp7",  file=FILE_SAL,  sheet="exp7 raw data",  type="temp_salinity", has_salinity=TRUE),
  list(id="exp8",  file=FILE_SAL,  sheet="exp8 raw data",  type="temp_salinity", has_salinity=TRUE),
  list(id="exp9",  file=FILE_SAL,  sheet="exp9 raw data",  type="temp_salinity", has_salinity=TRUE),
  list(id="exp10", file=FILE_SAL,  sheet="exp10 raw data", type="temp_salinity", has_salinity=TRUE),
  list(id="exp11", file=FILE_SAL,  sheet="exp11 raw data", type="temp_salinity", has_salinity=TRUE)
)

read_experiment <- function(exp) {
  df <- read_excel(exp$file, sheet=exp$sheet, .name_repair="universal") %>% clean_names()
  st_col  <- names(df)[str_detect(names(df), regex("surface|tension", ignore_case=TRUE))][1]
  od_col  <- names(df)[str_detect(names(df), regex("odaverage|od_average|od.average", ignore_case=TRUE))][1]
  sal_col <- names(df)[str_detect(names(df), regex("salinity", ignore_case=TRUE))][1]
  message("  ", exp$id, ": ST='", st_col, "' OD='", od_col, "'")
  df <- df %>% rename(st_raw=all_of(st_col), od_avg=all_of(od_col))
  if ("day" %in% names(df)) df <- df %>% rename(reading=day)
  df %>%
    mutate(
      experiment  = exp$id,
      exp_type    = exp$type,
      reading     = as.integer(reading),
      day         = day_map[as.character(reading)],
      strain      = as.character(strain),
      temperature = as.character(temperature),
      st_raw      = as.numeric(st_raw),
      od_avg      = as.numeric(od_avg),
      od1         = as.numeric(od1),
      od2         = as.numeric(od2),
      od3         = as.numeric(od3),
      salinity_m  = if (exp$has_salinity) as.numeric(.data[[sal_col]]) else 0
    ) %>%
    select(experiment, exp_type, strain, temperature, salinity_m,
           reading, day, st=st_raw, od=od_avg, od1, od2, od3)
}

message("Reading experiments...")
raw_df <- map(experiments, read_experiment) %>% bind_rows()
message("  Rows: ", nrow(raw_df), " | Experiments: ",
        n_distinct(raw_df$experiment), " | Strains: ", n_distinct(raw_df$strain))
message("  Temperature values found: ", paste(sort(unique(raw_df$temperature)), collapse=" | "))

message("\nQC checks...")
b <- raw_df %>% filter(day==1) %>%
  summarise(mn=mean(st,na.rm=T),lo=min(st,na.rm=T),hi=max(st,na.rm=T),sd=sd(st,na.rm=T))
message("  Baseline ST: mean=",round(b$mn,2)," range=[",round(b$lo,2),",",round(b$hi,2),"] SD=",round(b$sd,3))
if (b$lo<68.0||b$hi>71.0) warning("Baseline outside 68-71") else message("  Baseline: OK")
rc <- raw_df %>% group_by(experiment,strain,temperature,salinity_m) %>%
  summarise(n=n(),.groups="drop") %>% filter(n!=6)
if (nrow(rc)>0) {message("  WARNING unequal readings:"); print(rc)} else message("  Reading counts: OK")

message("\nComputing derived metrics...")
baseline_st <- raw_df %>% filter(day==1) %>%
  group_by(experiment,strain,temperature,salinity_m) %>%
  summarise(baseline_st=mean(st,na.rm=T),.groups="drop")

master <- raw_df %>%
  left_join(baseline_st, by=c("experiment","strain","temperature","salinity_m")) %>%
  mutate(
    delta_gamma    = baseline_st - st,
    delta_gamma_OD = if_else(od > OD_MIN_THRESHOLD, delta_gamma/od, NA_real_)
  )

message("\nComputing summary metrics...")
auc_trap <- function(days, values) {
  d <- data.frame(x=as.numeric(days), y=values) %>%
    filter(!is.na(y), !is.na(x)) %>% arrange(x)
  if (nrow(d)<2) return(NA_real_)
  sum(diff(d$x) * (head(d$y,-1) + tail(d$y,-1)) / 2)
}

summary_metrics <- master %>%
  group_by(experiment, exp_type, strain, temperature, salinity_m) %>%
  summarise(
    max_delta_gamma    = max(delta_gamma,    na.rm=TRUE),
    max_delta_gamma_OD = max(delta_gamma_OD, na.rm=TRUE),
    max_OD             = max(od,             na.rm=TRUE),
    final_st           = st[which.max(day)],
    final_delta_gamma  = delta_gamma[which.max(day)],
    AUC_delta_gamma    = auc_trap(day, delta_gamma),
    AUC_delta_gamma_OD = auc_trap(day, delta_gamma_OD),
    n_readings         = n(),
    .groups            = "drop"
  )
message("  Summary rows: ", nrow(summary_metrics))
message("  Temperature values in summary: ", paste(sort(unique(summary_metrics$temperature)), collapse=" | "))

message("\nAttaching strain metadata...")
master_full  <- master          %>% left_join(STRAIN_META, by=c("strain"="strain_phys"))
summary_full <- summary_metrics %>% left_join(STRAIN_META, by=c("strain"="strain_phys"))

no_wgs <- summary_full %>% filter(is.na(wgs_code)) %>% distinct(strain)
if (nrow(no_wgs)>0) message("  Physiology-only: ", paste(no_wgs$strain, collapse=", "))

message("\nComputing CSI...")
csi <- summary_full %>%
  filter(salinity_m==0, !is.na(AUC_delta_gamma_OD), !is.na(temperature)) %>%
  group_by(strain, temperature) %>%
  summarise(AUC_mean=mean(AUC_delta_gamma_OD, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=temperature, values_from=AUC_mean, values_fn=mean)

message("  CSI columns after pivot: ", paste(names(csi), collapse=" | "))

# Detect 4C and 30C columns by matching degree symbol variants
col4  <- names(csi)[str_detect(names(csi), "4")]
col30 <- names(csi)[str_detect(names(csi), "30")]
message("  4C col='", paste(col4,collapse=","), "' | 30C col='", paste(col30,collapse=","), "'")

if (length(col4)==1 && length(col30)==1) {
  csi <- csi %>%
    mutate(
      CSI       = .data[[col4]] / .data[[col30]],
      CSI_class = case_when(
        CSI > 1.5 ~ "Cold-stimulated",
        CSI < 0.7 ~ "Warm-preferred",
        TRUE      ~ "Neutral"
      )
    ) %>%
    select(strain, CSI, CSI_class)
  message("  Cold-stimulated: ", sum(csi$CSI_class=="Cold-stimulated",na.rm=TRUE))
  message("  Neutral:         ", sum(csi$CSI_class=="Neutral",na.rm=TRUE))
  message("  Warm-preferred:  ", sum(csi$CSI_class=="Warm-preferred",na.rm=TRUE))
  message("  NA (no 30C):     ", sum(is.na(csi$CSI)))
  print(csi %>% arrange(desc(CSI)))
} else {
  message("  WARNING: Could not detect 4C/30C columns. Showing all column names:")
  message("  ", paste(names(csi), collapse=" | "))
  csi <- csi %>% mutate(CSI=NA_real_, CSI_class=NA_character_) %>% select(strain,CSI,CSI_class)
}

summary_full <- summary_full %>% left_join(csi, by="strain")

message("\nSaving outputs...")
dir.create(DATA_PROC, recursive=TRUE, showWarnings=FALSE)
write_csv(master_full,  file.path(DATA_PROC, "master_dataset.csv"))
write_csv(summary_full, file.path(DATA_PROC, "summary_metrics.csv"))
write_csv(STRAIN_META %>% left_join(csi,by=c("strain_phys"="strain")), FILE_METADATA)
message("  master_dataset.csv:  ", nrow(master_full),  " rows")
message("  summary_metrics.csv: ", nrow(summary_full), " rows")
message("  strain_metadata.csv: written")

message("\n=== SANITY REPORT ===")
message("Total readings:       ", nrow(master_full))
message("Unique strains:       ", n_distinct(master_full$strain))
message("Experiments:          ", n_distinct(master_full$experiment))
message("delta_gamma range:    [",round(min(master_full$delta_gamma,na.rm=T),2),
        ", ",round(max(master_full$delta_gamma,na.rm=T),2),"] mN/m")
message("delta_gamma_OD range: [",round(min(master_full$delta_gamma_OD,na.rm=T),2),
        ", ",round(max(master_full$delta_gamma_OD,na.rm=T),2),"]")
message("\n=== MODULE 1 COMPLETE ===\n")
