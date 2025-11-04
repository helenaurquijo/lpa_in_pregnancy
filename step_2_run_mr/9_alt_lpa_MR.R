################################################################################
#          Run and plot MR analyses using alternative Lp(a) datasets           #
################################################################################

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

## Load env and paths
readRenviron("config.env")
projectdir <- Sys.getenv("projectdir")

## Get functions
source(paste0(projectdir, "/scripts/functions/MR_functions.R"))

################################################################################
# Prepare outcome data annotations

## Define outcome file path and list of outcomes of interest
out_path <- paste0(projectdir, "/data/outcomes")
out_info <- read_excel(paste0(projectdir,"/scripts/info.xlsx"), sheet = "out")

################################################################################
# Load exposure and outcome data

setwd(paste0(projectdir, "/results/sensitivity_analyses"))

## Define instrument file path
iv_path <- paste0(projectdir, "/data/exposures/instruments/")

## Load exposure files
exp_dat_mack <- fread(paste0(iv_path, "exp_dat_sens_lpa_mack2017.txt")) %>%
  as.data.frame() # already formated as exposure data

exp_dat_olink <- fread(paste0(iv_path, "exp_dat_sens_lpa_olink.txt")) %>%
  as.data.frame() # already formated as exposure data

## Combine the two exposure datasets
exp_dat <- bind_rows(exp_dat_mack, exp_dat_olink)

print(c("Exp file", dim(exp_dat))) 
# [1] "Exp file" "67"       "14"

## Load outcome file, keep relevant outcomes
out_dat <- 
  fread(paste0(out_path,"/main_out_dat.txt")) %>% 
  filter(outcome %in% out_info$outcome) %>%
  filter(SNP %in% exp_dat$SNP) %>%
  as.data.frame()

print(c("Out file", dim(out_dat)))
# "Out file" "930"     "16"

print(length(unique(out_dat$SNP)))
# [1] 65 ; 2 SNPs missing

## Harmonise, use eaf for palindromic SNPs
dat <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat,
  ) %>% 
  filter(mr_keep == TRUE) %>%
  mutate(
    study = "alt"
  )

print(c("Dat file", dim(dat)))
# [1] "Dat file" "949"      "36"

################################################################################
# Run main analysis

args <- dplyr::select(dat, study, outcome, exposure) %>%
  distinct %>%
  rename(stu = study, out = outcome, gene = exposure)

## Run MR
res <- 
  pmap(args, ~run_mr(..1, ..2, ..3, dat)) %>%
  bind_rows() %>%
  left_join(out_info %>% select(outcome, outcome_full, cat), by = "outcome") %>%
  mutate(
    # to reflect the effect of lowering Lp(a)
    b = -b
  ) %>%
  mutate(
    # to reflect the effect of lowering Lp(a) by 10nmol/L
    # 210nmol/L is approx 97.65mg/dL
    b_scaled = case_when(
      exposure == "LPA, apo(a) isoform-adjusted (mg/dL)" ~ b * 97.65,
      TRUE ~ b  # Default 
    ),
    se_scaled = case_when(
      exposure == "LPA, apo(a) isoform-adjusted (mg/dL)" ~ se * 97.65,
      TRUE ~ se  
    ),
  ) %>%
  calculate_or_ci() %>%
  mutate(
    FDR = p.adjust(pval, method = "fdr")
  )

## Save results
fwrite(
  res,
  file = "alt_lpa_mr_res.txt",
  sep = "\t",
  row.names = F, 
  col.names = T
  )

################################################################################
# Plot results

## Load plotting functions
source(paste0(projectdir, "/scripts/functions/plot_functions.R"))

## Load main results for comparison
main_res <- fread(paste0(projectdir, "/results/main_analysis/main_mr_res.txt")) %>%
  as.data.frame() %>%
  mutate(exposure = "LPA, immunoturbidimetric assay (main)")

to_plot <- bind_rows(main_res, res)

## Plot and save

colours <- qualitative_hcl(3, palette = "Dark 3")

suppressWarnings(
  plot_mr_special(to_plot, "exposure", scaled = TRUE) %>%
    ggsave(., 
      file = "alt_lpa_mr_res.png",
      width = 14,
      height = 8,
      dpi = 300)
)

q("no")