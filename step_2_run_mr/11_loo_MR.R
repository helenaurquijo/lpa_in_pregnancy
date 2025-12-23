################################################################################
#      Run and plot MR analyses adjusting leaving out one MR-PREG studies      #
################################################################################

##
# This script uses the study-specific GWAS estimates to generate study-specific
# MR estimates and meta-analyse them leaving out one study at a time.
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
library(metafor)
library(readxl)

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
# Load data

setwd(paste0(projectdir, "/results/sensitivity_analyses"))

## Load exposure data
iv_path <- paste0(projectdir,"/data/exposures/instruments/")

exp_dat <- fread(paste0(iv_path, "exp_dat_main_lpa.txt")) %>%
  as.data.frame() %>%
  format_data(
    type = "exposure",
    min_pval = .Machine$double.xmin
  ) %>%
  mutate(exposure = "Lp(a)")

## Load proxies info
proxies_info <- fread(paste0(iv_path, "proxies_info_main_lpa.txt"))

## Get outcome data by study and get proxies where needed
stu_outdat <- fread(paste0(out_path, "/stu_out_dat.txt")) %>% 
  dplyr::filter(Phenotype %in% out_info$outcome) %>%
  data.frame() %>% 
  split(., .$study) %>%
  map(., ~data.frame(.)) %>%
  map(., ~format_data(., type = "outcome")) %>%
  map(., ~mutate(., proxied = FALSE, proxy = NA))

stu_outdat <- stu_outdat %>%
  map(., function(out_dat) {
    map(unique(out_dat$outcome), ~replace_missing_with_proxies(out = .x, exp_dat, out_dat, proxies_info)) %>%
    bind_rows()
  })
  # N.B. conmal_any and pph not included as they only have one study

## Print outcomes for study
lapply(stu_outdat, function(x) {unique(x$outcome)})

################################################################################
# Run study-specific MR

## Harmonise data
studies <- names(stu_outdat) 
dat <- map(stu_outdat, ~harmonise_data(exp_dat, .)) %>%
  map2(., studies, ~mutate(.x, study = .y)) %>%
  bind_rows %>%
  filter(mr_keep == T)

## Set arguments for MR
args <- dplyr::select(dat, study, outcome, exposure) %>%
    distinct %>%
    rename(stu = study, out = outcome, gene = exposure)

## Run MR
stu_res <- 
  pmap(args, ~run_mr(..1, ..2, ..3, dat)) %>%
  bind_rows() %>%
  left_join(out_info %>% select(outcome, outcome_full, cat), by = "outcome") %>%
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
  stu_res,
  file = "stu_lpa_mr_res.txt",
  sep = "\t",
  row.names = F,
  col.names = T
)

################################################################################
# Run leave-one-out meta-analysis of MR results

## Set arguments for leave-one-out meta-analysis
loo_args <- expand_grid(
  gene = unique(stu_res$exposure),
  out = unique(stu_res$outcome)
)

## Run leave-one-out meta-analysis
loo_res <- pmap(loo_args, run_loo) %>%
  bind_rows %>%
  left_join(out_info %>% select(outcome, outcome_full, cat), by = "outcome") %>%
  arrange(match(outcome, out_info$outcome)) %>%
  rename(b = estimate) %>%
  mutate(
    # to reflect the effect of lowering Lp(a) by 210nmol/L
    b_scaled = b * 210,
    se_scaled = se * 210,
  ) %>%
  calculate_or_ci() %>%
  mutate(
    FDR = p.adjust(pval, method = "fdr")
  )

fwrite(
  loo_res,
  file = "loo_lpa_mr_res.txt",
  sep = "\t",
  row.names = F,
  col.names = T
)

################################################################################
# Plot leave-one-out results

## Load plotting functions
source(paste0(projectdir, "/scripts/functions/plot_functions.R"))

## Load main results for comparison
main_res <- fread(paste0(projectdir, "/results/main_analysis/main_mr_res.txt")) %>%
  mutate(study = "None") %>%
  # remove outcomes with one study only
  filter(!outcome %in% c("conmal_any", "pph")) %>% 
  as.data.frame()

## Merge data for plotting
to_plot <- loo_res %>%
  bind_rows(main_res) %>%
  arrange(match(outcome, out_info$outcome)) %>%
  rename(removing = study)

## Set plot colours  
colours <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#CC6677', '#882255', 'grey50')
names(colours) <- arrange(to_plot, removing) %>% pull(removing) %>% unique

## Plot and save
suppressWarnings(
  plot_mr_special(to_plot, "removing", scaled = TRUE) %>%
    ggsave(., 
      file = "loo_lpa_mr_res.png",
      width = 11,
      height = 11,
      dpi = 300)
)

q("no")