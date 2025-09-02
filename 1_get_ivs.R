################################################################################
#                Get LPA variant data for instrumental variables               #
################################################################################

##
# This scripts downloads, extracts and formats the 1) the original 29 SNPs from
# the Shi et al. 2024 paper which will be the main instruments 2) the SNPs from 
# the Mack et al. 2017 paper and 3) the SNPs from the Olink UKB-PPP GWAS. 
# SNPs from 2) and 3) will be used in one of the sensitivity analyses.
# Another script will identify proxies as not all of these are available in the 
# outcome datasets.
##

print(Sys.time())

library(data.table)
library(dplyr)
library(readxl)

## Load env and paths
readRenviron("scripts/config.env")
projectdir <- Sys.getenv("projectdir")

## Set up
setwd(paste0(projectdir, "/data/exposures/raw/"))

## 1) Shi et al. 2024

## Download the LPA variant data from medRxiv
# Evaluating genetically-predicted causal effects of lipoprotein(a) in human 
# diseases: a phenome-wide Mendelian randomization study. Shi et al. 2024,
# MedRxiv preprint.
system("wget https://www.medrxiv.org/content/medrxiv/early/2024/07/25/2024.07.24.24310950/DC3/embed/media-3.zip")

## Unzip the downloaded file
system("unzip -o media-3.zip")

## Extract summary statistics for finemapped LPA variants
dat <- read_excel("Supplementary_Tables/Supplementary_Table_2.xlsx")  %>%
  as.data.frame()

print(head(dat))

## Format the data for MR
exp_dat <- dat %>%
  filter(grepl("^rs", `rsID`)) %>% # filter out rows not start with rs
  rename(
    SNP = `rsID`,
    chr = `scaffold`,
    position = `position (GRCh38)`,
    other_allele = `non_effect_allele`,
    beta = `effect_size`,
    se = `stderr`,
    r2 = `variance_explained`,
    pval = `p-value`
  ) %>%
  mutate(chr = gsub("chr", "", chr)) #%>%
  #mutate(pval = pnorm(-abs(beta / se))) %>% # not sure if pointless as they still come out as zero
  #mutate(pval = pmax(pval, .Machine$double.xmin)) # is this necessary?

print(head(exp_dat))

## Save data
setwd(paste0(projectdir, "/data/exposures/instruments/"))

# 29 IVs from the paper
fwrite(exp_dat, file = "exp_dat_main_lpa.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# list of SNPs for proxy lookup
rsids <- data.frame(exp_dat$SNP)
fwrite(rsids, file = "rsids_main_lpa.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## 2) Mack et al. 2017

# (...)

## 3) Olink UKB-PPP GWAS

# (...)

q("no")
