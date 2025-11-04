################################################################################
#                    Genetic colocalization - Lp(a) x APPOs                    #
################################################################################

##
# This script runs the genetic colocalisation analyses assessing whether there
# is a shared causal cis variant for lp(a) levels and the APPOs at the LPA
# target locus
##

print(Sys.time())

## Initialise renv
source("renv/activate.R")

## Load packages
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(coloc)
library(readxl)

## Load env and paths
readRenviron("config.env")
projectdir <- Sys.getenv("projectdir")

## Get functions
source(paste0(projectdir, "/scripts/functions/coloc_functions.R"))

################################################################################
# Prepare outcome information

out_info <- read_excel(paste0(projectdir,"/scripts/info.xlsx"), sheet = "out")

################################################################################
# Prepare data


## Get lp(a) full summary stats
exp_file <- paste0(projectdir,"/data/exposures/full_sumstats/full_sumstats_lpa.txt")
exp_dat <- fread(exp_file) %>%
  # format for coloc to run well
  rename(MAF.exposure = eaf.exposure) %>%
  mutate(
    exposure = "Lp(a)",
    eaf.exposure = NA,
    samplesize.exposure = 253570
    )

print(dim(exp_dat)) # 3921 13

## Get and format APPO full summary stats

# From meta-analysis
out_file <- paste0(projectdir,"/data/outcomes/full_sumstats/ma_out_dat.txt")
outcome_dat_mrpeg <- as.data.frame(fread(out_file)) %>%
  filter(Phenotype %in% out_info$outcome) %>%
  format_data(type = "outcome") %>%
  split(.$outcome)

# PPH - Westergaard et al. 2024
out_file <- paste0(projectdir,"/data/outcomes/full_sumstats/pph_full_sumstats.txt")
outcome_dat_pph <- as.data.frame(fread(out_file)) %>%
  split(.$outcome)

# conmal_any - MoBa
out_file <- paste0(projectdir,"/data/outcomes/full_sumstats/conmal_full_sumstats.txt")
outcome_dat_conmal <- as.data.frame(fread(out_file)) %>%
  split(.$outcome)

## Combine all outcome datasets into one list
outcome_dat <- c(outcome_dat_mrpeg, outcome_dat_pph, outcome_dat_conmal)

lapply(outcome_dat, function(df) print(dim(df)))

## Create combinations of genes and tissues
combinations <- expand.grid(
  lipid = "Lp(a)",
  appo = out_info$outcome,
  stringsAsFactors = FALSE
)

################################################################################
# Run colocalization analyses

setwd(paste0(projectdir, "/results/sensitivity_analyses/"))

## Run coloc analysis 
coloc_results <- purrr::map2(
  combinations$lipid,
  combinations$appo,
  prep_and_run_coloc
) %>% list_rbind()

## See if any combinations did not run
combinations %>%
  mutate(
    trait1 = lipid,
    trait2 = appo) %>%
  anti_join(
    coloc_results %>%
    dplyr::select(trait1, trait2), by = c("trait1", "trait2"))
# [1] lipid  appo   trait1 trait2
# <0 rows> (or 0-length row.names)

## Save results
fwrite(
  coloc_results,
  "coloc.txt",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

################################################################################
# Plot local association plots

library(ieugwasr)
library(ggplot2)
library(geni.plots)

#source(paste0(projectdir, "/scripts/functions/MR_functions.R"))

## Format data for coloc and pass it to gassocplot function
dat_harmo <- purrr::map2(
  combinations$lipid,
  combinations$appo,
  format_coloc
)

## Create plots
walk(dat_harmo, ~run_assocplot(.))

q("no")