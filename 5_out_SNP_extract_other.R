################################################################################
#         Extract SNP-outcome summary data for for additional outcomes         #
################################################################################

##
# This scripts downloads (if required) and extracts SNP-outcome summary data for
# selected SNPs from the outcomes that are not part of the MR-PREG meta-analyses.
##

print(Sys.time())

library(data.table)
library(dplyr)
library(readxl)
library(TwoSampleMR)

## Load env and paths
readRenviron("scripts/config.env")
projectdir <- Sys.getenv("projectdir")

## Get outcome info
out_info <- read_excel(paste0(projectdir, "/scripts/info.xlsx"))

## Set up
setwd(paste0(projectdir, "/data/outcomes/raw"))

rsids <- fread(paste0(projectdir, "/scripts/rsids.txt"), header = FALSE)
nrow(rsids) #161 - includes proxies

## 1) PPH, from the GWAS meta-analysis of Westerg4ard et al. 2024 Nat Genet.

#     N.B. access needs to be obtained from 
#     https://download.decode.com/form/2024-steinthorsdottir-postpartu and 
#     the download cannot be performed programmatically

## Extract and format data
out_dat_pph <- fread("Postpartum_hemorrhage.txt.gz") %>%
  rename(
    SNP = rsName,
    chr = Chr,
    pos_b38 = PosB38,
    effect_allele = EA,
    other_allele = OA,
    eaf = EAfrq,
    beta = Effect,
    pval = P
  ) %>%
  filter(SNP %in% rsids$V1) %>%
  mutate(
    Phenotype = "pph",
    beta = log(beta), # effect size is given as OR
    se = abs(beta) / abs(qnorm(1 - pval / 2)), # se not provided
    ncase = 21512,
    ncontrol = 259500
    ) %>% 
  data.frame() %>%
  format_data(type = "outcome") %>%
  mutate(pos.outcome = NA)

print(head(out_dat_pph))

## 2) Congenital anomalies, from MoBa only so not part of the meta-analyses

## Read in, filter and format data
out_dat_conmal <- fread(Sys.getenv("conmal")) %>%
  mutate(Phenotype = "conmal_any") %>%
  filter(SNP %in% rsids$V1) %>%
  data.frame() %>%
  format_data(type = "outcome")

print(head(out_dat_conmal))

## Create a main out_dat file for all outcomes
ma_out_dat <- fread(paste0(projectdir, "/data/outcomes/ma_out_dat.txt")) %>%
  as.data.frame() %>%
  format_data(type = "outcome")

out_dat <- rbind(ma_out_dat, out_dat_pph, out_dat_conmal, fill = TRUE)

print(head(out_dat))

fwrite(out_dat, file = paste0(projectdir, "/data/outcomes/main_out_dat.txt"), sep = "\t")

q("no")