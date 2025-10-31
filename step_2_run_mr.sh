#!/bin/bash
#SBATCH --job-name=run_mr
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=3gb
#SBATCH --profile=task

# Exit immediately if any command fails
set -e

echo "Running script: $(basename $0)"
echo "Script run at: $(date '+%Y-%m-%d %H:%M:%S')"

source scripts/config.env

module load languages/R/4.3.3

## Run 8_main_MR.R
#  Runs main MR analysis - lp(a) x 21 APPOs
echo "Running 8_main_MR.R..."
R CMD BATCH ${projectdir}/scripts/step_2_run_mr/8_main_MR.R \
     ${projectdir}/logs/8_main_MR.Rout

## Run 9_alt_lpa_MR.R
#  Runs MR with alternative exposure datasets for lp(a) on the same outcomes
echo "Running 9_alt_lpa_MR.R..."
R CMD BATCH ${projectdir}/scripts/step_2_run_mr/9_alt_lpa_MR.R \
    ${projectdir}/logs/9_alt_lpa_MR.Rout

## Run 10_wlm_MR.R
#  Runs weighted linear model MR analysis
echo "Running 10_wlm_MR.R..."
R CMD BATCH ${projectdir}/scripts/step_2_run_mr/10_wlm_MR.R \
    ${projectdir}/logs/10_wlm_MR.Rout

## Run 11_loo_MR.R
#  Runs leave-one-out MR analysis for sensitivity testing
echo "Running 11_loo_MR.R..."
R CMD BATCH ${projectdir}/scripts/step_2_run_mr/11_loo_MR.R \
    ${projectdir}/logs/11_loo_MR.Rout

