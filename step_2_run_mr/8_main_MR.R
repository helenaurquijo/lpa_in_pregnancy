################################################################################
#                        Run and plot main MR analyses                         #
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

setwd(paste0(projectdir, "/results/main_analysis"))

## Define instrument file path
iv_path <- paste0(projectdir, "/data/exposures/instruments/")

## Load exposure files
exp_dat <- fread(paste0(iv_path, "exp_dat_main_lpa.txt")) %>%
  as.data.frame() %>%
  format_data(
    type = "exposure",
    min_pval = .Machine$double.xmin,
    # storing r2 as info column so it's not lost in format_data
    info_col = "r2" 
  ) %>%
  mutate(exposure = "Lp(a)")

print(c("Exp file", dim(exp_dat))) 
# [1] "Exp file" "29"       "13"

## Load outcome file, keep relevant outcomes
out_dat <- 
  fread(paste0(out_path,"/main_out_dat.txt")) %>% 
  filter(outcome %in% out_info$outcome) %>%
  mutate(
    proxied = FALSE,
    proxy = NA
  )

print(c("Out file", dim(out_dat)))
# [1] "Out file" "2938"     "18"

print(length(unique(out_dat$SNP)))
# [1] 219

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
# [1] "Dat file" "262"      "37"

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
  arrange(match(outcome, out_info$outcome)) %>%
  mutate(
    # to reflect the effect of lowering Lp(a)
    b = -b
  ) %>%
  mutate(
    # to reflect the effect of lowering Lp(a) by 210nmol/L
    b_scaled = b * 210,
    se_scaled = se * 210,
  ) %>%
  calculate_or_ci() %>%
  mutate(
    FDR = p.adjust(pval, method = "fdr")
  )

## Save results
fwrite(
  res,
  file = paste("main_mr_res.txt", sep = "_"),
  sep = "\t",
  row.names = F, 
  col.names = T
  )

################################################################################
# Plot results

## Load plotting functions
source(paste0(projectdir, "/scripts/functions/plot_functions.R"))

## Plot and save
suppressWarnings(
  plot_mr(res, scaled = TRUE) %>%
    ggsave(., 
      file = "main_mr_res_210nmol.png",
      width = 10,
      height = 8,
      dpi = 300)
)

suppressWarnings(
  plot_mr(res, scaled = FALSE) %>%
    ggsave(., 
      file = "main_mr_res_1nmol.png",
      width = 10,
      height = 8,
      dpi = 300)
)

q("no")