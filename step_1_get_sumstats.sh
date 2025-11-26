#!/bin/bash
#SBATCH --job-name=gets_ivs
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem=4gb
#SBATCH --profile=task

# Exit immediately if any command fails
set -e

echo "Running script: $(basename $0)"
echo "Script run at: $(date '+%Y-%m-%d %H:%M:%S')"

source scripts/config.env
cd scripts

module load languages/R/4.3.3

## Run 1_get_ivs.R
#  Extracts IV data from sources and formats
echo "Running 1_get_ivs.R..."
R CMD BATCH ${projectdir}/scripts/step_1_get_sumstats/1_get_ivs.R \
    ${projectdir}/logs/1_get_ivs.Rout

## Run 2_get_proxies.R
#  Searches proxies for main IV SNPs
echo "Running 2_get_proxies.R..."
R CMD BATCH ${projectdir}/scripts/step_1_get_sumstats/2_get_proxies.R \
    ${projectdir}/logs/2_get_proxies.Rout

## Run 3_out_SNP_extract.sh
#  Preps list of SNPs to extract from MR-PREG meta-analysis 
#  and calls the extraction R script (4_out_SNP_extract.R)
echo "Running 3_out_SNP_extract.sh..."
bash ${projectdir}/scripts/step_1_get_sumstats/3_out_SNP_extract.sh \
    > ${projectdir}/logs/3_out_SNP_extract.out 2>&1

## Run 5_out_SNP_extract_other.R
#  Extracts IV SNPs for outcomes that are not part of the MR-PREG analyses
echo "Running 5_out_SNP_extract_other.R..."
R CMD BATCH ${projectdir}/scripts/step_1_get_sumstats/5_out_SNP_extract_other.R \
    ${projectdir}/logs/5_out_SNP_extract_other.Rout

## Run 6_get_full_sumstats.R
#  Extracts full summary stats for Lp(a) for colocalization and MRlap analyses
echo "Running 6_get_full_sumstats.R..."
R CMD BATCH ${projectdir}/scripts/step_1_get_sumstats/6_get_full_sumstats.R \
    ${projectdir}/logs/6_get_full_sumstats.Rout

## Run 7_out_full_extract.sh
#  Extracts full summary statistics for outcome datasets
echo "Running 7_out_full_extract.sh..."
bash ${projectdir}/scripts/step_1_get_sumstats/7_out_full_extract.sh \
    > ${projectdir}/logs/7_out_full_extract.out 2>&1