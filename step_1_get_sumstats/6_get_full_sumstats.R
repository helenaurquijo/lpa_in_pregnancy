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

################################################################################
# Compare estimates to Shi et al. 2024 GWAS estimates (IV estimates)
# Sanity check to ensure estimates are highly correlated

library(ggplot2)
source(paste0(projectdir,"/scripts/functions/plot_functions.R"))

path <- paste0(projectdir, "/data/exposures/instruments/")
iv_dat <- fread(paste0(path, "exp_dat_main_lpa.txt"))

comp <- inner_join(formatted, iv_dat, by = "SNP") %>%
  select(
    SNP,
    effect_allele, effect_allele.exposure,
    other_allele, other_allele.exposure,
    eaf, MAF = eaf.exposure,
    beta, beta.exposure,
    se, se.exposure,
    pval, pval.exposure
    ) %>%
  mutate(
    # ensure no pvals of 0
    pval = ifelse(pval == 0, .Machine$double.xmin, pval),
    # Convert from SD units to nmol/L
    # SD as per UKB showcase https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30790
    beta.exposure = beta.exposure * 49.21
    )

## Check alleles are harmonised
allele_mismatch <- which(
  # Effect alleles do not match
  (comp$effect_allele != comp$effect_allele.exposure | 
   comp$other_allele != comp$other_allele.exposure) & 
  # And ensure they are merely flipped
  (comp$effect_allele == comp$other_allele.exposure | 
   comp$other_allele == comp$effect_allele.exposure)
)
if(length(allele_mismatch) == 0) {
  print("No allele mismatches found")
} else {
  print(paste("Found", length(allele_mismatch), "allele mismatches"))
}

## Ensure no palindromic SNPs
palindromic <- which((comp$effect_allele == "A" & comp$other_allele == "T") |
                     (comp$effect_allele == "T" & comp$other_allele == "A") |
                     (comp$effect_allele == "C" & comp$other_allele == "G") |
                     (comp$effect_allele == "G" & comp$other_allele == "C"))

palindromic
# the only palindromic SNP is harmonised and inferrable

## Flip alleles if needed
if(length(allele_mismatch) > 0) {
  # Fix each mismatched row directly
  for(i in allele_mismatch) {
    # Swap alleles
    comp$effect_allele[i] <- comp$effect_allele.exposure[i]
    comp$other_allele[i] <- comp$other_allele.exposure[i]
    # Flip the effect size
    comp$beta[i] <- -comp$beta[i]
  }
}

## Check alleles it has worked
allele_mismatch <- which(comp$effect_allele != comp$effect_allele.exposure | 
             comp$other_allele != comp$other_allele.exposure)
if(length(allele_mismatch) == 0) {
  print("No allele mismatches found")
} else {
  print(paste("Found", length(allele_mismatch), "allele mismatches"))
}

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