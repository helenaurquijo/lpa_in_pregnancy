#!/bin/bash

################################################################################
#               Extract SNP-outcome summary data for selected SNPs             #
################################################################################

##
# This scripts formats the final list of SNPs to be extracted fron the MR-PREG
# summmary statistics and calls the extraction R script. Outouts are moved to
# a specific folder so that when this script is called again it does not 
# overwrite this output.
##

echo "Running script: $(basename $0)"

module add languages/R/4.3.3

source scripts/config.env

tags=${projectdir}/data/exposures/instruments/main_lpa.tags
ivs=${projectdir}/data/exposures/instruments/rsids_main_lpa.txt

cat ${tags} ${ivs} | sort | uniq > ${projectdir}/scripts/rsids.txt

dos2unix ${projectdir}/scripts/rsids.txt

R CMD BATCH ${projectdir}/scripts/4_out_SNP_extract.R ${projectdir}/logs/4_out_SNP_extract.Rout 

# Move output files to avoid overwriting when extracting coloc SNPs
for f in ma_out_dat.txt stu_out_dat.txt duos_out_dat.txt trios_out_dat.txt; do
    mv "${projectdir}/data/$f" "${projectdir}/data/outcomes/$f"
done

# Remove temporary files
rm ${projectdir}/data/tmp_*