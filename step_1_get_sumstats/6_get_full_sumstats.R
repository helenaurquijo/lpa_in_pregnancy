################################################################################
#                     Extract full summary stats for Lp(a)                     #
################################################################################

##
# This script extracts full summary stats for coloc and MRlap analyses.
# These analyses require full summary data for all SNPs in the LPA locus and 
# those were not shared by Shi et al. 2024. Hence, we extract them from the
# GWAS by Sinnott-Armstrong et al. 2021 Nat Genet, also in the UK Biobank.
# 
# This script also extracts the list of SNPs for PPH and congenital anomalies
# outcomes, which are not part of the MR-PREG meta-analyses.
##

print(Sys.time())

## Initialise renv
source("renv/activate.R")

library(biomaRt)
library(data.table)
library(readxl)
library(TwoSampleMR)
library(Rsamtools)
library(dplyr)

## Load env and paths
readRenviron("config.env")
projectdir <- Sys.getenv("projectdir")

################################################################################
# Get gene window for LPA

setwd(paste0(projectdir, "/data/exposures/full_sumstats"))

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = "GRCh37"
  )

lpa_coords <- getBM(
  attributes = c("start_position", "end_position", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = "LPA",
  mart = mart
)

################################################################################
# Load full summary stats

#  Downloaded from https://nih.figshare.com/ndownloader/files/22769978
#    - Lipoprotein_A.imp.gz
#    - Lipoprotein_A.imp.gz.tbi
#  and saved in data/exposures/full_sumstats

## Query LPA region
file <- TabixFile("Lipoprotein_A.imp.gz")
start <- lpa_coords$start_position - 500000
end <- lpa_coords$end_position + 500000
query_res <- scanTabix(file, param = GRanges("6", IRanges(start, end))) %>% unlist()

## Convert to data frame
region_df <- fread(text = paste(query_res, collapse = "\n"))
colnames(region_df) <- c(
  "chr",
  "pos",
  "MarkerName",
  "other_allele",
  "effect_allele",
  "beta",
  "se",
  "pval",
  "Direction",
  "HetISq",
  "HetChiSq",
  "HetDf",
  "HetPVal",
  "MAF",
  "MACH_R2"
  )

print(nrow(region_df))

################################################################################
# Merge and format data

## Get UKB rsid metadata from chr6 
ukb_metadata <- Sys.getenv("ukb_metadata")
ukb_rsids <- fread(
  paste0(
    Sys.getenv("ukb_metadata"),
    "olink_rsid_map_mac5_info03_b0_7_chr6_patched_v2.tsv.gz")
  ) %>% 
  # select(1,4) %>%
  # Extract region of interest - GRCh37/hg19 build
  filter(POS19 >= start & POS19 <= end) %>%
  mutate(
    # Create new format: chr:pos19_ref_alt to match region_df
    MarkerName = paste0(6, ":", POS19, "_", REF, "_", ALT)
  ) %>% select(MarkerName, rsid)

## Merge with region data
merged <- left_join(region_df, ukb_rsids, by = "MarkerName") %>%
  select(SNP = rsid, chr, pos, other_allele, effect_allele, beta, se, pval, MAF) %>%
  as.data.frame()

print(nrow(merged))

## Format data
formatted <- format_data(
  merged,
  type = "exposure",
  eaf_col = "MAF",
  min_pval = .Machine$double.xmin
) %>%
  mutate(exposure = "Lp(a)")

## Save data
fwrite(
  formatted,
  "full_sumstats_lpa.txt",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
  )

## Save list of SNPs and set up for the SNP-out extraction script
fwrite(
  data.frame(SNP = formatted$SNP),
  "full_rsids.txt",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE
  )

################################################################################
# Compare estimates to Shi et al. 2024 GWAS estimates (IV estimates)
# Sanity check to ensure estimates are highly correlated

library(ggplot2)
source(paste0(projectdir,"/scripts/functions/plot_functions.R"))

## Extract main IVs and format as outcome for harmonization
path <- paste0(projectdir, "/data/exposures/instruments/")
iv_dat <- fread(paste0(path, "exp_dat_main_lpa.txt")) %>%
  as.data.frame() %>%
  format_data(., type = "outcome")

## Harmonise data
comp <- harmonise_data(
    exposure_dat = formatted,
    outcome_dat = iv_dat,
    action = 2
  ) %>%
  mutate(
    # ensure no pvals of 0
    pval.outcome = ifelse(pval.outcome == 0, .Machine$double.xmin, pval.outcome),
    pval.exposure = ifelse(pval.exposure == 0, .Machine$double.xmin, pval.exposure),
    # Convert from SD units to nmol/L
    # SD as per UKB showcase https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30790
    beta.exposure = beta.exposure * 49.21
  )

## Plot comparison
png("lpa_gwas_vs_iv.png", 
        res = 400, width = 25, height = 25, units = "cm")
scatter_plot(comp)
dev.off()


################################################################################
# Extract LPA region for outcomes not in MR-PREG meta-analyses

setwd(paste0(projectdir, "/data/outcomes/full_sumstats"))

## 1) PPH, from the GWAS meta-analysis of Westerg4ard et al. 2024 Nat Genet.
pph_path <- paste0(projectdir, "/data/outcomes/raw/")
out_dat_pph <- fread(paste0(pph_path,"Postpartum_hemorrhage.txt.gz")) %>%
  rename(
    SNP = rsName,
    chr = Chr,
    pos_b38 = PosB38,
    effect_allele = EA,
    other_allele = OA,
    eaf = EAfrq,
    OR = Effect,
    pval = P
  ) %>%
  # filter for SNPs in Lp(a) GWAS as you need data for the same SNPs
  filter(SNP %in% formatted$SNP) %>%
  mutate(
    Phenotype = "pph",
    beta = log(OR), # effect size is given as OR
    se = abs(beta) / abs(qnorm(1 - pval / 2)), # se not provided
    ncase = 21512,
    ncontrol = 259500
    ) %>% 
  data.frame() %>%
  format_data(type = "outcome") %>%
  mutate(pos.outcome = NA)

print(dim(out_dat_pph))

## Save data
fwrite(
  out_dat_pph,
  "pph_full_sumstats.txt",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
  )

## 2) Congenital anomalies, from MoBa only so not part of the meta-analyses

## Read in, filter and format data
out_dat_conmal <- fread(Sys.getenv("conmal")) %>%
  mutate(Phenotype = "conmal_any") %>%
  # filter for SNPs in Lp(a) GWAS as you need data for the same SNPs
  filter(SNP %in% formatted$SNP) %>%
  data.frame() %>%
  format_data(type = "outcome")

print(dim(out_dat_conmal))

## Save data
fwrite(
  out_dat_conmal,
  "conmal_full_sumstats.txt",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
  )

q("no")