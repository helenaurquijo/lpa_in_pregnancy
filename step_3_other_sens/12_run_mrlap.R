################################################################################
#             Run and plot MR analyses adjusting for sample overlap            #
################################################################################

##
# This script runs MR using MRlap to adjust for sample overlap between the
# exposure dataset (entirely UKB) and the outcome datasets (some partly UKB).
# TO DO: think of way to make this more efficient as it repeats the LDSC for
# LP(a) every time, script is requiring almost 50GB of RAM to run.
##

print(Sys.time())

## Initialise renv
source("renv/activate.R")

## Load packages
library(TwoSampleMR)
library(MendelianRandomization)
library(ieugwasr)
library(tidyverse)
library(data.table)
library(ggforestplot)
library(gridExtra)
library(patchwork)
library(readxl)
library(colorspace)
library(MRlap)

## Load env and paths
readRenviron("config.env")
projectdir <- Sys.getenv("projectdir")
mrpregdir <- Sys.getenv("mrpregdir")

## Get functions
source(paste0(projectdir, "/scripts/functions/MR_functions.R"))

## Set wd
setwd(paste0(projectdir, "/results/sensitivity_analyses"))

################################################################################
# Prepare outcome data annotations

## Define outcome file path and list of outcomes of interest
out_path <- paste0(projectdir, "/data/outcomes")
out_info <- read_excel(paste0(projectdir,"/scripts/info.xlsx"), sheet = "out")

################################################################################
# Prepare downloaded LD and hm3 files for MRlap

# Downloaded from https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v

## Unzip folder - run only once
system(
  paste0("unzip ", projectdir, "/data/for_mrlap/eur_w_ld_chr.zip ",
  "-d ", projectdir, "/data/for_mrlap/")
  )

################################################################################
# Load MR-PREG liftover output
#  - output from liftover performed for 21099247 variants that were present in 
#    the maternal and/or offspring MR PREG meta-analysed summary statistics. 
#    Since these include UKB and FinnGen this cover most of the variants needed.

liftover_dat <- fread(
  paste0(mrpregdir,"/results/liftover/mrpreg_rsids_22m_dbSnp155_hg19_hg38.txt.gz"),
  sep = " "
)

################################################################################
# Load exposure data and format

## Define lp(a) full summary stats path
lpa_file <- paste0(projectdir,"/data/exposures/full_sumstats/Lipoprotein_A.imp.gz")

# Read in full GWAS summary stats and merge with liftover output
exp_dat <- fread(lpa_file) %>%
  rename(
    chr = "#CHROM",
    pos_hg19 = POS,
  ) %>%
  mutate(
    chr = as.character(chr),
    N = 253570
  ) %>%
  left_join(., liftover_dat, by = c("chr", "pos_hg19"), relationship = "many-to-many") %>%
  mutate(Phenotype = "Lp(a)") %>%
  rename(
    CHR = chr,
    POS = pos_hg19,
    SNP = rsid,
    BETA = Effect,
    SE = StdErr,
  ) %>% 
  data.frame()

print(dim(exp_dat)) # 8888991 19

################################################################################
# Get and format APPO full summary stats

# From meta-analysis
meta_info <- read_csv(paste0(projectdir,"/scripts/vars_metanalysis_v5_mum.csv"))
mrpreg_outs <- filter(meta_info, Variable_name %in% out_info$outcome)
mrpreg_full <- list()
for (out in mrpreg_outs$Variable_name) {
  print(out)
  df <- fread(paste0(Sys.getenv("mrpregrawgwas"), "metanalyses-R3.mum.", out, ".txt.gz")) %>%
    rename(
      CHR = chr,
      POS = pos,
      ALT = effect_allele,
      REF = other_allele,
      EAF = eaf,
      BETA = beta,
      SE = se,
      P = pval,
      N = samplesize
    ) %>%
    mutate(
      outcome = out
    ) %>%
    data.frame()
  mrpreg_full[[out]] <- df
}

# PPH - Westergaard et al. 2024 - needs grch37 coordinates
pph_file <- paste0(projectdir, "/data/outcomes/raw/Postpartum_hemorrhage.txt.gz")

out_dat_pph <- fread(pph_file) %>%
  mutate(Chr = gsub("chr", "", Chr)) %>% # Convert from "chr1" to "1"
  rename(
    rsid = rsName,
    chr = Chr,
    pos_hg38 = PosB38,
    ALT = EA,
    REF = OA,
    EAF = EAfrq,
    OR = Effect,
  ) %>%
  left_join(., liftover_dat, by = c("rsid", "chr", "pos_hg38"), relationship = "many-to-many") %>%
  mutate(
    Phenotype = "pph",
    BETA = log(OR), # effect size is given as OR
    SE = abs(BETA) / abs(qnorm(1 - P / 2)), # se not provided
    N = 281012
  ) %>% 
  rename(
    CHR = chr,
    POS = pos_hg19,
    SNP = rsid,
  ) %>%
  data.frame()

# conmal_any - MoBa
out_file <- paste0(projectdir,"/data/outcomes/full_sumstats/conmal_full_sumstats.txt")
out_dat_conmal <- fread(Sys.getenv("conmal")) %>%
  rename(
    CHR = chr,
    POS = pos,
    ALT = effect_allele,
    REF = other_allele,
    EAF = eaf,
    BETA = beta,
    SE = se,
    P = pval,
    N = samplesize
  ) %>%
  mutate(Phenotype = "conmal_any") %>%
  data.frame()

## Combine all outcome datasets into one list
out_dat <- mrpreg_full 
out_dat[['pph']] <- out_dat_pph 
out_dat[['conmal_any']] <- out_dat_conmal

lapply(out_dat, function(df) print(dim(df)))

## Set function arguments
args <- expand.grid(
  exp = unique(exp_dat$Phenotype),
  out = names(out_dat)
) 

################################################################################
# Run MRlap

res <- pmap(args, function(exp, out) {
  
  message(paste0("Running MRlap for ", exp, " and ", out))
  
  # Subset exposure and outcome data
  exp_sub <- subset(exp_dat, Phenotype == exp)
  out_sub <- out_dat[[out]]

  # Calculate MAF only if it doesn't already exist
  if (!"MAF" %in% colnames(exp_sub)) {
    exp_sub$MAF <- pmin(exp_sub$EAF, 1 - exp_sub$EAF)
  }

  if (!"MAF" %in% colnames(out_sub)) {
    out_sub$MAF <- pmin(out_sub$EAF, 1 - out_sub$EAF)
  }
  
  # Calculate Z if not already present
  if (!"Z" %in% colnames(exp_sub)) {
    exp_sub$Z <- exp_sub$BETA / exp_sub$SE
  }
  if (!"Z" %in% colnames(out_sub)) {
    out_sub$Z <- out_sub$BETA / out_sub$SE
  }
  
  # Run MRlap
  mrlap_res <- tryCatch({
    MRlap(
      exposure = exp_sub,
      outcome = out_sub,
      # downloaded from https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v
      ld = paste0(projectdir, "/data/for_mrlap/eur_w_ld_chr/"),
      hm3 = paste0(projectdir, "/data/for_mrlap/eur_w_ld_chr/w_hm3.snplist"),
      save_logfiles = TRUE,
      verbose = TRUE
    )
  }, error = function(e) {
    warning(paste("MRlap failed for", exp, "and", out, ":", e$message))
    return(NULL)
  })
  
  if (is.null(mrlap_res)) {
    return(NULL)
  }
  
  # Extract results from MRlap output
  res_df <- data.frame(
    exposure = exp,
    outcome = out,
    obs_b = mrlap_res$MRcorrection$observed_effect,
    obs_se = mrlap_res$MRcorrection$observed_effect_se,
    obs_p = mrlap_res$MRcorrection$observed_effect_p,
    corrected_b = mrlap_res$MRcorrection$corrected_effect,
    corrected_se = mrlap_res$MRcorrection$corrected_effect_se,
    corrected_p = mrlap_res$MRcorrection$corrected_effect_p,
    test_difference = mrlap_res$MRcorrection$test_difference,
    p_difference = mrlap_res$MRcorrection$p_difference,
    n_IVs = mrlap_res$MRcorrection$m_IVs,
    IVs = paste(mrlap_res$MRcorrection$IVs, collapse = "|")  # Join IVs with semicolon
  )
  
  return(res_df)

}) %>%
  bind_rows()

## Format results
formatted_res <- res %>%
  # Add any additional formatting here if needed
  mutate(
    # Convert effects to reflect lowering Lp(a)
    obs_b = -obs_b,
    corrected_b = -corrected_b
  ) %>%
  mutate(
    # to reflect the effect of lowering Lp(a) by 210nmol/L
    obs_b_scaled = obs_b * 120,
    corrected_b_scaled = corrected_b * 120,
    obs_se_scaled = obs_se * 120,
    corrected_se_scaled = corrected_se * 120
  )

## Save results
fwrite(
  formatted_res,
  file = "mrlap_results.txt",
  sep = "\t",
  row.names = FALSE, 
  col.names = TRUE
)