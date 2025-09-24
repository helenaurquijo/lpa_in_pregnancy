################################################################################
#          Run and plot MR analyses using alternative Lp(a) datasets           #
################################################################################

rm(list = ls())

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
readRenviron("scripts/config.env")
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
exp_dat <- rbind(exp_dat_mack, exp_dat_olink)

print(c("Exp file", dim(exp_dat))) 
# [1] "Exp file" "29"       "12"

## Load outcome file, keep relevant outcomes
out_dat <- 
  fread(paste0(out_path,"/main_out_dat.txt")) %>% 
  filter(outcome %in% out_info$outcome)

print(c("Out file", dim(out_dat)))
# 

print(length(unique(out_dat$SNP)))
# 

## Harmonise, use eaf for palindromic SNPs
dat <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat,
  ) %>% 
  filter(mr_keep == TRUE) %>%
  mutate(
    study = "main"
  )

print(c("Dat file", dim(dat)))
# [1] "Dat file" "245"      "36"

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
    # to reflect the effect of lowering Lp(a) by 10nmol/L
    # scaling needs to be discussed
    b = -b,
    se = se,
  ) %>%
  mutate(
    # Conditionally calculate OR, UCI, and LCI for binary outcomes
    OR = if_else(cat == "bin", exp(b), NA_real_),
    UCI = if_else(cat == "bin", exp(b + 1.96 * se), NA_real_),
    LCI = if_else(cat == "bin", exp(b - 1.96 * se), NA_real_),
    FDR = p.adjust(pval, method = "fdr")
  )

## Save results
fwrite(
  res,
  file = paste("alt_lpa_mr_res.txt", sep = "_"),
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
  rename(exposure = "Lp(a), immunoturbidimetric assay (main)")

to_plot <- rbind(main_res, res)

## Plot and save

plot_mr(res) %>%
  ggsave(., 
    file = "alt_lpa_mr_res.png",
    width = 12,
    height = 8,
    dpi = 300)

q("no")