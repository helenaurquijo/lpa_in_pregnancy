#!/bin/bash
#SBATCH --job-name=other_sens
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=40gb
#SBATCH --profile=task

# Exit immediately if any command fails
set -e

echo "Running script: $(basename $0)"
echo "Script run at: $(date '+%Y-%m-%d %H:%M:%S')"

source config.env

module load languages/R/4.3.3

## Run 12_run_mrlap.R
#  Analyzes sample overlap between exposure and outcome datasets
echo "Running 12_run_mrlap.R..."
R CMD BATCH ${projectdir}/scripts/step_3_other_sens/12_run_mrlap.R \
    ${projectdir}/logs/12_run_mrlap.Rout

## Run 13_coloc.R
#  Performs colocalization analysis
echo "Running 13_coloc.R..."
R CMD BATCH ${projectdir}/scripts/step_3_other_sens/13_coloc.R \
    ${projectdir}/logs/13_coloc.Rout