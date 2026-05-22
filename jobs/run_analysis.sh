#!/bin/bash
#SBATCH --job-name=physiology_analysis
#SBATCH --output=../results/logs/%j_analysis.log
#SBATCH --error=../results/logs/%j_analysis.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# Load R module (adjust version to what is available on IOCAS HPC)
module load R/4.3.0

cd /data/habib/Physiology

# Run modules in sequence
Rscript scripts/01_build_master_dataset.R
Rscript scripts/02_exploratory_analysis.R
Rscript scripts/03_temperature_effects.R
Rscript scripts/04_salinity_temperature_interaction.R
Rscript scripts/05_strain_classification.R
Rscript scripts/06_genomic_integration.R
Rscript scripts/07_manuscript_figures.R
