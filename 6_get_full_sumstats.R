################################################################################
#                     Extract full summary stats for Lp(a)                     #
################################################################################

##
# This script extracts full summary stats for coloc and MRlap analyses.
# These analyses require full summary data for all SNPs in the LPA locus and 
# those were not shared by Shi et al. 2024. Hence, we extract them from the
# GWAS by Sinnott-Armstrong et al. 2021 Nat Genet, also in the UK Biobank.
# This script also extracts the list of SNPs for PPH and congenital anomalies
# outcomes, which are not part of the MR-PREG meta-analyses.
##

print(Sys.time())

library(biomaRt)
library(data.table)
library(readxl)
library(TwoSampleMR)
library(Rsamtools)
library(dplyr)

## Load env and paths
readRenviron("scripts/config.env")
projectdir <- Sys.getenv("projectdir")

################################################################################
# Get gene window for LPA

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

setwd(paste0(projectdir, "/data/exposures/full_sumstats"))

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

## Get rsids
ukb_metadata <- Sys.getenv("ukb_metadata")
ukb_rsids <- fread(
  paste0(
    Sys.getenv("ukb_metadata"),
    "olink_rsid_map_mac5_info03_b0_7_chr6_patched_v2.tsv.gz")
  ) %>% 
  select(1,4) %>%
  mutate(
    # Extract components from 6:63979:C:T:imp:v1 format in first column
    chr = gsub(":.*", "", .[[1]]),
    pos = gsub(".*?:(\\d+):.*", "\\1", .[[1]]),
    ref = gsub(".*?:\\d+:(\\w+):(\\w+):.*", "\\1", .[[1]]),
    alt = gsub(".*?:\\d+:(\\w+):(\\w+):.*", "\\2", .[[1]]),
    # Create new format: chr:pos_ref_alt
    MarkerName = paste0(chr, ":", pos, "_", ref, "_", alt)
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

q("no")