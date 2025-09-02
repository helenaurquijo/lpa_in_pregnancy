#!/bin/bash

################################################################################
#                          Get proxies for Lp(a) IVs                           #
################################################################################

##
# This script aims to extract a list of tagging variants for each of the 
# finamapped variants from the original finemapping Lp(a) SNPs paper.
##

module load plink/1.9-beta6.27-openblas
source scripts/config.env

echo "Running script: $(basename $0)"
echo "Script run at: $(date '+%Y-%m-%d %H:%M:%S')"

refpan=${mrpregdir}/results/UKB_10k/ukb_eur_10k
rsids=${projectdir}/data/exposures/instruments/rsids_main_lpa.txt
window=500
r2=0.8

plink --bfile ${refpan} \
  --show-tags ${rsids} \
  --list-all \
  --tag-kb ${window} \
  --tag-r2 ${r2} \
  --out ${projectdir}/data/exposures/instruments/main_lpa

