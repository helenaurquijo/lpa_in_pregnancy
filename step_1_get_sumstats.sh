#!/bin/bash
#SBATCH --job-name=gets_ivs
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=200mb
#SBATCH --profile=task

echo "Running script: $(basename $0)"
echo "Script run at: $(date '+%Y-%m-%d %H:%M:%S')"

source scripts/config.env

module load languages/R/4.3.3

## Run 1_get_ivs.R
#  Extracts IV data from sources and formats

R CMD BATCH ${projectdir}/scripts/1_get_ivs.R ${projectdir}/logs/1_get_ivs.Rout

## Run 2_get_proxies.R
#  Searches proxies for main IV SNPs 

R CMD BATCH ${projectdir}/scripts/2_get_proxies.R > ${projectdir}/logs/2_get_proxies.Rout

## Run 3_out_SNP_extract.sh
#  Preps list of SNPs to extract from MR-PREG meta-analysis 
#  and calls the extraction R script (#4)

bash ${projectdir}/scripts/3_out_SNP_extract.sh > ${projectdir}/logs/3_out_SNP_extract.out 

## Run 5_out_SNP_extract_other.R
#  Extracts IV SNPs for outcomes that are not part pf the MR-PREG analyses

R CMD BATCH ${projectdir}/scripts/5_out_SNP_extract_other.R ${projectdir}/logs/5_out_SNP_extract_other.Rout