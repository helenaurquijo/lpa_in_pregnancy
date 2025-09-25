#!/bin/bash
#SBATCH --job-name=run_mr
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=200mb
#SBATCH --profile=task

# Exit immediately if any command fails
set -e

echo "Running script: $(basename $0)"
echo "Script run at: $(date '+%Y-%m-%d %H:%M:%S')"

source scripts/config.env

module load languages/R/4.3.3

## Run 6_main_MR.R
#  Runs main MR analysis - lp(a) x 21 APPOs
echo "Running 6_main_MR.R..."
R CMD BATCH ${projectdir}/scripts/6_main_MR.R ${projectdir}/logs/6_main_MR.Rout

## Run 7_alt_lpa_MR.R
#  Runs MR with alternative exposure datasets for lp(a) on the same outcomes
echo "Running 7_alt_lpa_MR.R..."
R CMD BATCH ${projectdir}/scripts/7_alt_lpa_MR.R ${projectdir}/logs/7_alt_lpa_MR.Rout