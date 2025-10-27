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

library(biomaRt)
library(data.table)
library(dplyr)
library(readxl)
library(readr)
library(vroom)
library(TwoSampleMR)

## Load env and paths
readRenviron("scripts/config.env")
projectdir <- Sys.getenv("projectdir")

## Set up
setwd(paste0(projectdir, "/data/exposures/raw/"))

## Get gene window for LPA
# GRCh38 coordinates from Ensembl Biomart
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = "115" # Set 2025 version - https://sep2025.archive.ensembl.org
  )

lpa_coords <- getBM(
  attributes = c("start_position", "end_position", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = "LPA",
  mart = mart
)

###############################################################################
# 1) Shi et al. 2024

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
  mutate(`rsID` = gsub(":.*", "", `rsID`)) %>% # to improve matching with ref panel 
  rename(
    SNP = `rsID`,
    chr = `scaffold`,
    position = `position (GRCh38)`,
    other_allele = `non_effect_allele`,
    beta = `effect_size`,
    se = `stderr`,
    eaf = `effect_allele_freq`,
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

###############################################################################
# 2) Mack et al. 2017

setwd(paste0(projectdir, "/data/exposures/raw/"))

# The summary statistics from this GWAS were shared in PDF format
# https://www.jlr.org/cms/10.1194/jlr.M076232/attachment/cfbb7466-d752-4ab7-b82b-558b638259e2/mmc1.pdf
# The information had to be manually extracted from ST8 and columns renamed as below:

# SNP: SNP
# effect_allele: Min_allele
# other_allele:	Ref_allele
# direction: Effect Direction	
# From the Single SNP analysis, fixed effects model:
# s_beta_int: beta(se) INT scale
# s_beta: beta(se) original scale
# s_pval:	p-value
# From the Joint SNP analysis, fixed effects model, adjusted for all 31 independent SNPs:
# j_beta_int:	beta(se) INT scale
# j_beta:	beta(se) original scale
# j_pval: p-value

mack_dat <- read_csv("mack_2017_lpa_snps.csv", ) %>%
  as.data.frame() %>%
  mutate(
    s_se_int = gsub(".*\\(", "", s_beta_int),
    s_se_int = gsub("\\)", "", s_se_int),
    j_se_int = gsub(".*\\(", "", j_beta_int),
    j_se_int = gsub("\\)", "", j_se_int),
    s_se = gsub(".*\\(", "", s_beta),
    s_se = gsub("\\)", "", s_se),
    j_se = gsub(".*\\(", "", j_beta),
    j_se = gsub("\\)", "", j_se),
    s_beta_int = gsub("\\(.*", "", s_beta_int),
    j_beta_int = gsub("\\(.*", "", j_beta_int),
    s_beta = gsub("\\(.*", "", s_beta),
    j_beta = gsub("\\(.*", "", j_beta),
  ) %>%
  # force numeric
  mutate_at(
    vars(
      s_beta_int, s_se_int, s_beta, s_se,
      j_beta_int, j_se_int, j_beta, j_se,
      s_pval, j_pval
    ), as.numeric
  )

# EAF (MAF as effect allele is the minor allele) was also provided in ST4 and 
# columns were renamed as below:
# SNP: SNP
# `chr:pos`: chr:position
# near_gene: nearest gene
# effect_allele: Minor allele
# other_allele:	ref allele
# maf_cohorts: MAF range of contributing studies
# maf_1kg3: MAF 1000G ph1v3

mack_info <- fread("mack_2017_snp_info.csv") %>%
  as.data.frame() %>%
  select(SNP, maf_1kg3, `chr:pos`) %>%
  # separate `chr:pos` into chr and position
  tidyr::separate(`chr:pos`, into = c("chr", "pos"), sep = ":") %>%
  # force numeric
  mutate_at(vars(chr, pos), as.numeric) %>%
  rename(eaf = maf_1kg3)

## Merge the two datasets and format for MR
#  Using marginal estimates (as opposed to joint model) to match the outcome 
#  GWAS estimates, in the original scale (as opposed to IRNT) to match 
#  the approach with the main exposure dataset better

exp_dat_mack <- merge(mack_dat, mack_info, by = "SNP", all.x = TRUE) %>%
  rename(
    beta = s_beta,
    se = s_se,
    pval = s_pval
  ) %>%
  select(SNP, effect_allele, other_allele, eaf, beta, se, pval, chr, pos) %>%
  mutate(Phenotype = "LPA, apo(a) isoform-adjusted (mg/dL)") %>%
  format_data(.,
    type = "exposure",
    min_pval = .Machine$double.xmin    
  ) 

print(head(exp_dat_mack))
nrow(exp_dat_mack) # [1] 30

## Save data
setwd(paste0(projectdir, "/data/exposures/instruments/"))

fwrite(exp_dat_mack, file = "exp_dat_sens_lpa_mack2017.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# List of SNPs for extraction from outcome datasets
rsids_mack <- data.frame(exp_dat_mack$SNP)
fwrite(rsids_mack, file = "rsids_sens_lpa_mack2017.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Save conditional estimates (joint model) for use as weights in the 
# genetic score analysis
f_name <- paste0(projectdir, "/data/UKB_RAP/lpa_mack2017_weights.txt")
weights <- merge(mack_dat, mack_info, by = "SNP", all.x = TRUE) %>%
  rename(
    beta = j_beta,
    se = j_se,
    pval = j_pval
  ) %>%
  select(SNP, effect_allele, other_allele, eaf, beta, se, pval, chr, pos)
fwrite(weights, file = f_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################
# 3) Olink UKB-PPP GWAS
#    Summary statistics can be downloaded from:
#    https://www.synapse.org/Synapse:syn51364943

## Set up to extract the relevant .gz file from the tar file
setwd(paste0(projectdir, "/data/exposures/raw/"))
dir_name <- "LPA_P08519_OID30747_v1_Inflammation_II"
file_name <- "discovery_chr6_LPA:P08519:OID30747:v1:Inflammation_II.gz"
dir.create("tmp", showWarnings = FALSE)

# Extract file from tar
tar_file <- file.path(Sys.getenv("ukb_pqtl_dir"), paste0(dir_name, ".tar"))
untar(tar_file, files = file.path(dir_name, file_name), exdir = "tmp")
  
# Move and cleanup
file.rename(
  file.path("tmp", dir_name, file_name),
  file.path("tmp", file_name)
)

unlink(file.path("tmp", dir_name), recursive = TRUE)

## Read in the data
lpa_olink <- vroom(file.path("tmp", file_name), col_select = -c(7, 9, 12, 14))
  
colnames(lpa_olink) <- c("chr", "pos_hg38", "ID", "other_allele", "effect_allele", "eaf", "samplesize", "beta", "se", "log_pval")

## Data does not have rsids - load the mapping file to get rsids
#  This is found in the metadata folder in the same link as above
key_p <- vroom(
    paste0(
      Sys.getenv("ukb_metadata"),
      "olink_rsid_map_mac5_info03_b0_7_chr6_patched_v2.tsv.gz"), 
    col_select = c(1,4)
    )

# Merge to get rsids, filter and format
lpa_olink_dat <- merge(lpa_olink, key_p, by = "ID") %>% 
  mutate(Phenotype = "LPA, Olink UKB-PPP") %>%
  rename(SNP = rsid) %>%
  select(-ID) %>%
  # filter to p < 5e-8
  filter(log_pval > -log10(5e-8)) %>%
  # restrict to SNPs within vicinity of the LPA gene
  filter(
    pos_hg38 >= (lpa_coords$start_position - 500000) & pos_hg38 <= (lpa_coords$start_position + 500000)
    ) %>%
  # format for MR
  format_data(.,
    type = "exposure",
    chr_col = "chr",
    pos_col = "pos_hg38",
    min_pval = .Machine$double.xmin    
  )

print(head(lpa_olink_dat))
nrow(lpa_olink_dat) # 3321

## Clump to get indendent SNPs (r2 < 0.01, 1000kb)
#  Need to discuss the parameters here
lpa_olink_exp_dat <- clump_data(
  lpa_olink_dat,
  clump_r2 = 0.01,
  clump_kb = 1000,
  bfile = Sys.getenv("refpan"),
  plink_bin = plinkbinr::get_plink_exe()
  )

dim(lpa_olink_exp_dat) # 37 SNPs

## Save data
setwd(paste0(projectdir, "/data/exposures/instruments/"))

fwrite(lpa_olink_exp_dat, file = "exp_dat_sens_lpa_olink.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# List of SNPs for extraction from outcome datasets
rsids_olink <- data.frame(lpa_olink_exp_dat$SNP)
fwrite(rsids_olink, file = "rsids_sens_lpa_olink.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Final cleanup - remove tmp directory
setwd(paste0(projectdir, "/data/exposures/raw/"))
unlink("tmp", recursive = TRUE)

q("no")
