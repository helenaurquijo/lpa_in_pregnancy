#!/bin/bash

module add languages/R/4.3.3

source scripts/config.env

tags=${projectdir}/data/exposures/instruments/main_lpa.tags
ivs=${projectdir}/data/exposures/instruments/rsids_main_lpa.txt

cat ${tags} ${ivs} | sort | uniq > ${projectdir}/scripts/rsids.txt

dos2unix ${projectdir}/scripts/rsids.txt

R CMD BATCH ${projectdir}/scripts/3_b_out_SNP_extract.R ${projectdir}/logs/3_b_out_SNP_extract.Rout 

# Move output files to avoid overwriting when extracting coloc SNPs
for f in ma_out_dat.txt stu_out_dat.txt duos_out_dat.txt trios_out_dat.txt; do
    mv "${projectdir}/data/$f" "${projectdir}/data/outcomes/$f"
done

# Remove temporary files
rm ${projectdir}/data/tmp_*
rm ${projectdir}/scripts/tmp_*
