# Biosurfactant Physiology Analysis
### Temperature and salinity effects on biosurfactant production in Qaidam Basin and marine isolates

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R version](https://img.shields.io/badge/R-%3E%3D4.3.0-blue)](https://www.r-project.org/)

---

## Overview

This repository contains the complete statistical analysis pipeline for Chapter 3 of the PhD thesis:

> **"Biosurfactant production under combined temperature and salinity stress in Qaidam Basin and marine bacterial isolates: a genotype-phenotype analysis"**

The study investigates whether low temperature enhances biosurfactant production in isolates from the Qaidam Basin (Tibet, China), and whether this effect is modulated by high salinity — testing the hypothesis that cold-adapted biosurfactant production is a specialized ecological adaptation linked to specific biosynthetic gene clusters (BGCs).

**Related repository:** [WGS-Biosurfactant-Pipeline](https://github.com/Habib-Ullah1/WGS-Biosurfactant-Pipeline) — whole genome sequencing and BGC mining pipeline for the same isolates (cited in Module 6).

---

## Experimental design

| Factor | Levels |
|---|---|
| Temperature | 4°C · 15°C · 30°C |
| Salinity | 0.5M · 1.0M · 1.5M NaCl (salinity experiments only) |
| Strains | 22 unique isolates (17 with WGS, 5 physiology-only) |
| Origin | 15 Qaidam Basin · 7 marine |
| Time points | 6 readings per experiment (Days 1, 3, 5, 7, 9, 11) |
| Experiments | 10 total (exp2–exp6: temperature-only · exp7–exp11: temperature × salinity) |

**Response variables:** Surface tension (mN/m), OD (growth), Δγ (ST reduction), Δγ/OD (biomass-normalised production efficiency), AUC(Δγ/OD) (sustained production over 11 days).

---

## Repository structure

```
.
├── scripts/
│   ├── 00_config.R                         # paths, constants, colour palettes
│   ├── 01_build_master_dataset.R           # data harmonisation, Δγ/OD, AUC
│   ├── 02_exploratory_analysis.R           # descriptive stats, producer classification
│   ├── 03_temperature_effects.R            # LMM: temperature main effect (H1)
│   ├── 04_salinity_temperature_interaction.R  # LMM: Temp × Salinity (H2)
│   ├── 05_strain_classification.R          # clustering, OD vs production
│   ├── 06_genomic_integration.R            # genotype-phenotype linkage (H3)
│   └── 07_manuscript_figures.R             # final publication-quality figures
├── data/
│   ├── processed/
│   │   ├── master_dataset.csv              # main tidy dataset (all experiments)
│   │   └── strain_metadata.csv            # strain × WGS × species × BGC table
│   └── wgs_inputs/                         # TSVs from WGS pipeline (symlinks on HPC)
│       ├── bgc_type_matrix.tsv
│       ├── stress_genes_all.tsv
│       ├── cog_matrix.tsv
│       └── strain_metadata_all21.tsv
├── results/
│   ├── module02_exploratory/
│   ├── module03_temperature/
│   ├── module04_interaction/
│   ├── module05_classification/
│   └── module06_genomic/
├── figures/
│   ├── main/                               # Fig 1–5 (manuscript)
│   └── supplementary/                      # Supplementary figures
├── tables/
│   ├── main/                               # Tables 1–3 (manuscript)
│   └── supplementary/
├── docs/
│   └── statistical_methods.md             # detailed methods for manuscript
└── jobs/
    └── run_analysis.sh                     # SLURM job script (IOCAS HPC)
```

---

## Analysis modules

| Module | Script | Description | Key output |
|---|---|---|---|
| 1 | `01_build_master_dataset.R` | Harmonise all 10 experiments, compute Δγ/OD, AUC | `master_dataset.csv` |
| 2 | `02_exploratory_analysis.R` | Descriptive stats, Gaussian mixture producer classification | Classification table, baseline check |
| 3 | `03_temperature_effects.R` | LMM: Temperature ~ Δγ/OD + (1\|Strain) + (1\|Experiment) | F-stat, pairwise, Cohen's d |
| 4 | `04_salinity_temperature_interaction.R` | LMM: Temp × Salinity interaction, Bliss Independence | Interaction plots, synergy/antagonism |
| 5 | `05_strain_classification.R` | Hierarchical clustering, cold stimulation index, OD scatter | Heatmap, cluster assignments |
| 6 | `06_genomic_integration.R` | BGC presence × Max(Δγ/OD) correlation, stress gene linkage | Genotype-phenotype figures |
| 7 | `07_manuscript_figures.R` | Assemble all publication-ready multi-panel figures | Fig 1–5 + Supplementary |

---

## Key hypotheses

- **H1 (Cold stimulation):** Mean Δγ/OD at 4°C > 15°C > 30°C across strains (main effect of temperature)
- **H2 (Salinity modulation):** Cold-stimulation effect is amplified at higher salinity (Temperature × Salinity interaction)
- **H3 (Strain heterogeneity):** Strains with lipopeptide BGCs (Lichenysin/Surfactin) show significantly higher Max(Δγ/OD) at 4°C vs strains without

---

## Requirements

### R packages
```r
# Install all dependencies
install.packages(c(
  # Data wrangling
  "tidyverse", "readxl", "janitor",
  # Statistical modelling
  "lme4", "lmerTest", "emmeans", "MuMIn", "mclust", "coin",
  # Visualisation
  "ggplot2", "patchwork", "cowplot", "ggbeeswarm",
  "pheatmap", "ComplexHeatmap", "RColorBrewer", "viridis",
  # Reporting
  "broom.mixed", "knitr", "rmarkdown"
))
```

### R version
Developed and tested with R ≥ 4.3.0

---

## Reproducibility

All scripts source `scripts/00_config.R` at the top, which sets all file paths and random seeds. To reproduce the full analysis:

```bash
# On IOCAS HPC
cd /data/habib/Physiology
Rscript scripts/01_build_master_dataset.R
Rscript scripts/02_exploratory_analysis.R
# ... or submit all at once:
sbatch jobs/run_analysis.sh
```

To reproduce locally (off-HPC), update the `WGS_RESULTS` path in `scripts/00_config.R` to point to local copies of the WGS downstream files.

---

## Data availability

Raw Excel files (`Temp_experiments20260520.xlsx`, `Salinity_and_Temp_Combined_experiments_data.xlsx`) are not tracked in this repository (unpublished experimental data). The derived `data/processed/master_dataset.csv` is included and sufficient to reproduce all statistical results and figures.

WGS-derived input files (BGC matrix, stress genes, COG profiles) are outputs of the companion pipeline at [WGS-Biosurfactant-Pipeline](https://github.com/Habib-Ullah1/WGS-Biosurfactant-Pipeline).

---

## Citation

> [To be updated upon manuscript acceptance]

---

## Author

**Habib Ullah**
Institute of Oceanology, Chinese Academy of Sciences (IOCAS)
PhD Candidate, Microbial Ecology

