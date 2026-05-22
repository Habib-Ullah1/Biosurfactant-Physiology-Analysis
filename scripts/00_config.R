# =============================================================================
# 00_config.R
# Central configuration — sourced at the top of every analysis script.
# Edit ONLY this file when paths change (e.g., running off-HPC).
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Reproducibility
# -----------------------------------------------------------------------------
set.seed(2024)

# -----------------------------------------------------------------------------
# 2. Paths — HPC
# Change WGS_RESULTS to a local path if running off-HPC
# -----------------------------------------------------------------------------
PROJECT_ROOT  <- "/data/habib/Physiology"
DATA_RAW      <- file.path(PROJECT_ROOT, "data/raw")
DATA_PROC     <- file.path(PROJECT_ROOT, "data/processed")
WGS_INPUTS    <- file.path(PROJECT_ROOT, "data/wgs_inputs")
WGS_RESULTS   <- "/data/habib/WGS/analysis/11_downstream"   # companion pipeline

RESULTS       <- file.path(PROJECT_ROOT, "results")
FIGS_MAIN     <- file.path(PROJECT_ROOT, "figures/main")
FIGS_SUPP     <- file.path(PROJECT_ROOT, "figures/supplementary")
TABLES_MAIN   <- file.path(PROJECT_ROOT, "tables/main")
TABLES_SUPP   <- file.path(PROJECT_ROOT, "tables/supplementary")

# Raw data files
FILE_TEMP <- file.path(DATA_RAW, "Temp_experiments20260520.xlsx")
FILE_SAL  <- file.path(DATA_RAW, "Salinity_and_Temp_Combined_experiments_data.xlsx")

# Processed master dataset (written by Module 1, read by all others)
FILE_MASTER   <- file.path(DATA_PROC, "master_dataset.csv")
FILE_METADATA <- file.path(DATA_PROC, "strain_metadata.csv")

# WGS-derived inputs for Module 6
FILE_BGC_MATRIX   <- file.path(WGS_INPUTS, "bgc_type_matrix.tsv")
FILE_BGC_COMPLETE <- file.path(WGS_INPUTS, "bgc_complete_table.tsv")
FILE_STRESS_GENES <- file.path(WGS_INPUTS, "stress_genes_all.tsv")
FILE_COG_MATRIX   <- file.path(WGS_INPUTS, "cog_matrix.tsv")
FILE_WGS_META     <- file.path(WGS_INPUTS, "strain_metadata_all21.tsv")

# -----------------------------------------------------------------------------
# 3. Experimental constants
# -----------------------------------------------------------------------------

# Mapping: Reading index (1-6) -> actual day
DAY_MAP <- c("1" = 1, "2" = 3, "3" = 5, "4" = 7, "5" = 9, "6" = 11)

# Temperatures (ordered factor for plots)
TEMP_LEVELS  <- c("4C", "15C", "30C")

# Salinities (ordered factor for plots)
SAL_LEVELS   <- c(0.5, 1.0, 1.5)

# Minimum OD to compute Δγ/OD (below this = unreliable normalisation)
OD_MIN_THRESHOLD <- 0.01

# Minimum Δγ to classify as any production (noise floor)
DELTA_GAMMA_MIN  <- 2.0   # mN/m

# Baseline surface tension expected range (for QC check)
ST_BASELINE_EXPECTED <- c(68.5, 71.0)  # mN/m

# -----------------------------------------------------------------------------
# 4. Strain metadata — complete cross-reference table
# WGS code | Physiology strain name | Species | Environment | Lipopeptide BGC
# -----------------------------------------------------------------------------
STRAIN_META <- data.frame(
  strain_phys = c(
    "SW1C-HB20", "SW2C3",   "SW4C6",   "SW1B1",    "SW2B2",
    "H2E-10C",   "H5G",     "4-2C3",   "M2-3-i9",  "H5B",
    "SW2B-HB1",  "SW2BHB1", "4-2A1",   "C8-P4",
    "SX203",     "KR141",   "KR108",
    "4-2C6",     "SW1B2",   "SW1C2",   "M2-6B2",   "SX206"
  ),
  wgs_code = c(
    "A1",  "A2",  "A3",  "A4",  "A5",
    "A6",  "A7",  "A8",  "A9",  "A10",
    "A11", "A12", "A13", "A14",
    "L1",  "L2",  "L3",
    NA,    NA,    NA,    NA,    NA
  ),
  species = c(
    "B. zhangzhouensis",    "B. zhangzhouensis",    "B. zhangzhouensis",
    "B. zhangzhouensis",    "B. zhangzhouensis",    "P. frigoritolerans",
    "B. thuringiensis",     "B. zhangzhouensis",    "B. zhangzhouensis",
    "B. thuringiensis",     "Peribacillus sp.",      "B. pumilus",
    "Virgibacillus sp.",    "B. zhangzhouensis",
    "Sphingobium sp.",      "Sphingobium sp.",       "Kocuria rosea",
    NA, NA, NA, NA, NA
  ),
  environment = c(
    rep("Qaidam", 14),
    rep("Marine", 3),
    rep(NA, 5)
  ),
  has_lipopeptide_bgc = c(
    TRUE,  TRUE,  TRUE,  TRUE,  TRUE,   # A1-A5: Lichenysin/Surfactin
    FALSE, FALSE, TRUE,  TRUE,  FALSE,  # A6-A10
    FALSE, TRUE,  FALSE, TRUE,          # A11-A14
    FALSE, FALSE, FALSE,                # L1-L3: no lipopeptide BGC
    NA, NA, NA, NA, NA                  # no WGS
  ),
  bgc_type_lipopeptide = c(
    "Lichenysin", "Surfactin",  "Lichenysin", "Lichenysin", "Lichenysin",
    NA,           NA,           "Lichenysin", "Lichenysin", NA,
    NA,           "Lichenysin", NA,           "Lichenysin",
    NA,           NA,           NA,
    NA, NA, NA, NA, NA
  ),
  has_wgs = c(
    rep(TRUE, 17),
    rep(FALSE, 5)
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 5. Colour palettes — consistent across ALL figures
# -----------------------------------------------------------------------------

# Temperature colours (blue = cold, grey = ambient, red = warm)
COL_TEMP <- c(
  "4C"  = "#378ADD",   # blue
  "15C" = "#888780",   # grey
  "30C" = "#D85A30"    # coral/red
)

# Salinity colours (light to dark teal)
COL_SAL <- c(
  "0.5" = "#9FE1CB",
  "1"   = "#1D9E75",
  "1.5" = "#085041"
)

# Environment colours
COL_ENV <- c(
  "Qaidam" = "#534AB7",   # purple
  "Marine" = "#378ADD"    # blue
)

# Producer classification colours
COL_PRODUCER <- c(
  "High"    = "#1D9E75",   # teal
  "Moderate"= "#EF9F27",   # amber
  "Low"     = "#D3D1C7"    # light grey
)

# Cold stimulation classification colours
COL_CSI <- c(
  "Cold-stimulated" = "#378ADD",
  "Neutral"         = "#888780",
  "Warm-preferred"  = "#D85A30"
)

# -----------------------------------------------------------------------------
# 6. ggplot2 theme — publication-ready base theme
# Applied to every figure via: + THEME_PUB
# -----------------------------------------------------------------------------
library(ggplot2)

THEME_PUB <- theme_classic(base_size = 11) +
  theme(
    # Text
    axis.title       = element_text(size = 11, colour = "black"),
    axis.text        = element_text(size = 10, colour = "black"),
    strip.text       = element_text(size = 10, face = "bold"),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    # Panels
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.5),
    # Legend
    legend.background = element_blank(),
    legend.key        = element_blank(),
    # Margins
    plot.margin = margin(8, 8, 8, 8, "pt")
  )

# Figure dimensions for Microbiome journal
# Single column: 88mm wide; Double column: 180mm wide; Max height: 230mm
FIG_W_SINGLE <- 88  / 25.4  # inches
FIG_W_DOUBLE <- 180 / 25.4  # inches
FIG_H_MAX    <- 230 / 25.4  # inches

# Resolution for saving
FIG_DPI <- 300

# -----------------------------------------------------------------------------
# 7. Helper function — save figures consistently
# Usage: save_fig("module03_temperature/fig_temp_main", plot_obj, w=FIG_W_DOUBLE, h=4)
# -----------------------------------------------------------------------------
save_fig <- function(name, plot, w = FIG_W_DOUBLE, h = 4,
                     dir = FIGS_MAIN, dpi = FIG_DPI) {
  path_pdf <- file.path(dir, paste0(name, ".pdf"))
  path_png <- file.path(dir, paste0(name, ".png"))
  dir.create(dirname(path_pdf), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path_pdf, plot, width = w, height = h, dpi = dpi)
  ggplot2::ggsave(path_png, plot, width = w, height = h, dpi = dpi)
  message("Saved: ", path_pdf)
  message("Saved: ", path_png)
}

# -----------------------------------------------------------------------------
# 8. Helper function — save tables consistently
# -----------------------------------------------------------------------------
save_table <- function(df, name, dir = TABLES_MAIN) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".csv"))
  write.csv(df, path, row.names = FALSE)
  message("Saved: ", path)
}

# -----------------------------------------------------------------------------
# 9. Startup message
# -----------------------------------------------------------------------------
message("=== 00_config.R loaded ===")
message("Project root : ", PROJECT_ROOT)
message("R version    : ", R.version$version.string)
message("Seed         : 2024")

