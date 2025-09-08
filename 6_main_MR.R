################################################################################
#                        Run and plot main MR analyses                         #
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
# Run main analysis

setwd(paste0(projectdir, "/results/main_analysis"))

## Define instrument file path
iv_path <- paste0(projectdir, "/data/exposures/instruments/")

## Load exposure files
exp_dat <- fread(paste0(iv_path, "exp_dat_main_lpa.txt")) %>%
  as.data.frame() %>%
  format_data(
    type = "exposure",
    min_pval = .Machine$double.xmin
  ) %>%
  mutate(exposure = "LPA")

print(c("Exp file", dim(exp_dat))) 
# [1] "Exp file" "29"       "12"

## Load outcome file, keep relevant outcomes
out_dat <- 
  fread(paste0(out_path,"/main_out_dat.txt")) %>% 
  filter(outcome %in% out_info$outcome) %>%
  mutate(
    proxied = FALSE,
    proxy = NA
  )

print(c("Out file", dim(out_dat)))
# [1] "Out file" "1899"     "18"

print(length(unique(out_dat$SNP)))
# [1] 162

## Find proxies where needed
proxies_info <- fread(paste0(iv_path, "proxies_info_main_lpa.txt"))

out_dat <- map(unique(out_dat$outcome), ~replace_missing_with_proxies(out = .x, exp_dat, out_dat, proxies_info)) %>%
  bind_rows()

table(out_dat$proxied)

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

args <- dplyr::select(dat, study, outcome, exposure) %>%
  distinct %>%
  rename(stu = study, out = outcome, gene = exposure)

## order the outcomes 
order <- out_info$outcome

## Run MR
res <- 
  pmap(args, ~run_mr(..1, ..2, ..3, dat)) %>%
  bind_rows() %>%
  left_join(out_info %>% select(outcome, cat), by = "outcome") %>%
  mutate(
    b = -b, # to reflect the effect of lowering Lp(a)
    # Conditionally calculate OR, UCI, and LCI for binary outcomes
    OR = if_else(cat == "bin", exp(b), NA_real_),
    UCI = if_else(cat == "bin", exp(b + 1.96 * se), NA_real_),
    LCI = if_else(cat == "bin", exp(b - 1.96 * se), NA_real_),
    outcome = factor(outcome, levels = order)
  )

