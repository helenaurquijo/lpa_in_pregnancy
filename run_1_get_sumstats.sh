#!/bin/bash
#SBATCH --job-name=gets_ivs
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --mem=5gb
#SBATCH --profile=task

echo "Script run at: $(date '+%Y-%m-%d %H:%M:%S')"

source scripts/config.env

load languages/R/4.3.3

## Run 1_a_get_ivs.R
#  Extracts IV data from sources and formats

R CMD BATCH ${projectdir}/scripts/1_a_get_ivs.R ${projectdir}/logs/1_a_get_ivs.Rout

## Run 1_b_get_proxies.sh
#  Searches proxies for main IV SNPs 

bash ${projectdir}/scripts/1_b_get_proxies.sh ${projectdir}/logs/1_b_get_proxies.out

## Run 3_a_out_SNP_extract.sh
#  Preps list of SNPs to extract from MR-PREG meta-analysis 
#  and calls the extraction R script

bash ${projectdir}/scripts/3_a_out_SNP_extract.sh ${projectdir}/logs/3_a_out_SNP_extract.out