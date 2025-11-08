#!/bin/bash

################################################################################
#                Extract SNP-outcome summary data for LPA locus                #
################################################################################

##
# This scripts formats the full list of SNPs in the LPA region to be extracted 
# from the MR-PREG summmary statistics and calls the extraction R script. 
# Outputs should not overwrite the output from when the script was run earlier
# as they were moved elsewhere.
##

# Exit immediately if any command fails
set -e

echo "Running script: $(basename $0)"

module add languages/R/4.3.3

source config.env

snps=${projectdir}/data/exposures/full_sumstats/full_rsids.txt

cat ${snps} | sort | uniq > ${projectdir}/scripts/rsids.txt

dos2unix ${projectdir}/scripts/rsids.txt

R CMD BATCH ${projectdir}/scripts/step_1_get_sumstats/4_out_SNP_extract.R ${projectdir}/logs/4_full_out_SNP_extract.Rout 

# Move output files to avoid overwriting when extracting coloc SNPs
for f in ma_out_dat.txt stu_out_dat.txt duos_out_dat.txt trios_out_dat.txt; do
    mv "${projectdir}/data/$f" "${projectdir}/data/outcomes/full_sumstats/$f"
done

# Remove temporary files
rm ${projectdir}/data/tmp_*