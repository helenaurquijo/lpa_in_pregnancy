################################################################################
#                  Run and plot positive control MR analyses                   #
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

setwd(paste0(projectdir, "/results/iv_checks"))

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
# [1] "Exp file" "29"       "12"

## Extract CAD outcome data from OpenGWAS
out_dat <-
  extract_outcome_data(
    snps = exp_dat$SNP,
    outcomes = "ieu-a-7",
    proxies = 1, 
    opengwas_jwt = get_opengwas_jwt()
  ) %>% mutate(outcome = "CAD")

dim(out_dat)
#[1] 13 16

## Harmonise, use eaf for palindromic SNPs
dat <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat,
  ) %>% 
  filter(mr_keep == TRUE) %>%
  mutate(study = "iv_check")

print(c("Dat file", dim(dat)))
# [1] "Dat file" "13"      "35"

################################################################################
# Run analysis

res <- 
  run_mr("iv_check", "CAD", "Lp(a)", dat) %>%
  mutate(cat = "bin") %>%
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
  file = paste("pos_ctrl_mr_res.txt", sep = "_"),
  sep = "\t",
  row.names = F, 
  col.names = T
  )

print(res)

################################################################################
# Plot, comparing to the main analysis results

source(paste0(projectdir, "/scripts/functions/plot_functions.R"))

res1 <- res %>% 
  mutate(
    outcome_full = "Coronary artery disease",
    cat = "bin"
  )

main_res <- fread(paste0(projectdir, "/results/main_analysis/main_mr_res.txt"))

to_plot <- bind_rows(res1, main_res)

suppressWarnings(
  plot_mr(to_plot, scaled = TRUE) %>%
    ggsave(., 
      file = "pos_ctrl_mr_res_210nmol.png",
      width = 10,
      height = 8,
      dpi = 300)
)

suppressWarnings(
  plot_mr(to_plot, scaled = FALSE) %>%
    ggsave(., 
      file = "pos_ctrl_mr_res_1nmol.png",
      width = 10,
      height = 8,
      dpi = 300)
)

q("no")