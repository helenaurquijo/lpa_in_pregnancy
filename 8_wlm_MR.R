################################################################################
#         Run and plot MR analyses adjusting for fetal genetic effects         #
################################################################################

##
# This script uses the wighted linear model (WLM) conditional GWAS estimates 
# and performs MR analyses adjusting for effects on fetal genotype on the
# outcomes, which may constitute horizontal pleiotropy.
##

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
# Load data

setwd(paste0(projectdir, "results/sensitivity_analyses"))

## Load exposure data
iv_path <- paste0(projectdir,"data/exposures/instruments/")

exp_dat <- fread(paste0(iv_path, "exp_dat_main_lpa.txt")) %>%
  as.data.frame() %>%
  format_data(
    type = "exposure",
    min_pval = .Machine$double.xmin
  ) %>%
  mutate(exposure = "Lp(a)")

## Load outcome file
duo_dat <- fread(paste0(out_path, "/duos_out_dat.txt")) %>% 
  filter(!is.na(se_mat_donuts)) # Duos
trio_dat <- fread(paste0(out_path, "/trios_out_dat.txt")) %>%
  filter(!is.na(se_mat_donuts)) # Trios

## Load WLM conmal_any data
conmal_wlm <- Sys.getenv("conmalwlm")
rsids <- fread(paste0(iv_path, "rsids.txt"), header = FALSE)

conmal_duo_dat <- fread(paste0(conmal_wlm, "conmal_any_donuts_duos_mum-child.txt.gz")) %>%
  mutate(
    Phenotype = "conmal_any",
    type = "duos_mum-child",
    analyses = "duos_mum-child",
    study = "MOBA") %>%
  filter(SNP %in% rsids$V1) %>%
  data.frame()

conmal_trio_dat <- fread(paste0(conmal_wlm, "conmal_any_donuts_trios.txt.gz")) %>%
  mutate(
    Phenotype = "conmal_any",
    type = "trios",
    analyses = "trios",
    study = "MOBA") %>%
  filter(SNP %in% rsids$V1) %>%
  data.frame()

## Append conmal to duo and trio datasets
duo_dat <- rbind(duo_dat, conmal_duo_dat)
trio_dat <- rbind(trio_dat, conmal_trio_dat)

## Format outcome data
format_outdat <- function(df, type) {
  dat <- data.frame(df) %>%
  format_data(
    type = "outcome",
    snps = exp_dat$SNP,
    eaf_col = "eaf_mat",
    beta_col = paste0("beta_", type),
    se_col = paste0("se_", type),
    pval_col = paste0("p_", type)
    ) %>%
    filter(outcome %in% out_info$outcome)
}

m <- c(
  "Maternal (unadjusted)", "Maternal (fetal adjusted)", "Maternal (fetal+paternal adjusted)",
  "Fetal (unadjusted)", "Fetal (maternal adjusted)", "Fetal (maternal+paternal adjusted)",
  "Paternal (unadjusted)", "Paternal (maternal+fetal adjusted)"
)

t <- c(
  "mat", "mat_donuts", "mat_donuts",
  "off", "off_donuts", "off_donuts",
  "pat", "pat_donuts"
  )

d <- list(
  duo_dat, duo_dat, trio_dat,
  duo_dat, duo_dat, trio_dat,
  trio_dat, trio_dat
  )

## Format and name outcome data
out_dat <- map2(d, t, ~format_outdat(df = .x, type = .y)) 
names(out_dat) <- m

## Harmonise data 
dat <- map(out_dat, ~harmonise_data(exp_dat, .)) %>%
  map2(., names(out_dat), ~mutate(.x, adjustment = .y)) %>%
  bind_rows()

################################################################################
# Run MR analysis

## Create list of args for the MR run
args <- select(dat, exposure, outcome, adjustment) %>%
  distinct %>%
  rename(gene = exposure, out = outcome, adj = adjustment)

## Run MR with WLM estimates
wlm_res <- pmap(args, function( out, adj) {

  print(paste("Lp(a)", "-", out, "-", adj))

  ## Subset data
  df <- filter(
    dat,
    outcome == out &
    adjustment == adj &
    mr_keep == TRUE
  )
  print(head(df))

  ## empty df if N SNP = 0
  if(nrow(df)==0) {

    print("No SNPs selected")

    res2 <- data.frame(
      outcome = out,
      exposure = unique(df$exposure),
			adjustment = adj,
      nsnp = 0
		)

  } else if (nrow(df)==1) {

    ## Run MR Wald ratio if N SNP = 1
    res <- mr(df)
    print("Running Wald ratio")

    res2 <- mutate(res,
      exposure = unique(df$exposure),
      adjustment = adj,
      method = "Wald ratio"
    )

  } else if(nrow(df)>1) {

    ## Run MR IVW accounting for LD if N SNPs > 1
    print("Running IVW")
    df2 <- dat_to_MRInput_me(df, get_correlation = TRUE)

    if(!is.null(df2[[1]])) {

      res <- MendelianRandomization::mr_ivw(df2[[1]])     
      res2 <- data.frame(
        exposure = unique(df$exposure),
        adjustment = adj,
        exposure = res$Exposure,
        outcome = res$Outcome,
        b = res$Estimate,
        se = res$StdError,
        pval = res$Pvalue,
        nsnp = res$SNPs,
        het_Q = res$Heter.Stat[1],
        het_pval = res$Heter.Stat[2],
        method = "Generalised IVW"
      )

      print("gIVW successful")

    } else {

      print("Could not be harmonised - LD clump and run IVW or Wald ratio")

      res <- df %>%
        mutate(rsid = SNP, pval = pval.exposure) %>%
        ld_clump_local(.,
          clump_kb = 10000,
          clump_r2 = 0.01,
          clump_p = 0.99,
          bfile = paste0(mrpregdir, "results/UKB_10k/ukb_eur_10k"),
          plink_bin = plinkbinr::get_plink_exe()
        ) %>%
        mr

      res2 <- res %>%
        mutate(
          exposure = unique(df$exposure),
          adjustment = adj,
          method = ifelse(nsnp == 1, "Wald ratio", "IVW")
        )
    }
  }

  res2$study <- "wlm"
  return(res2)
  
}) %>%
  bind_rows %>%
  mutate(
    # to reflect the effect of lowering Lp(a) by 10nmol/L
    b = -b * 10,
    se = se * 10,
  ) %>%
  mutate(
    # Conditionally calculate OR, UCI, and LCI for binary outcomes
    OR = if_else(cat == "bin", exp(b), NA_real_),
    UCI = if_else(cat == "bin", exp(b + 1.96 * se), NA_real_),
    LCI = if_else(cat == "bin", exp(b - 1.96 * se), NA_real_),
    FDR = p.adjust(pval, method = "fdr")
  )

fwrite(
  wlm_res,
  file = "wlm_lpa_mr_res.txt",
  sep = "\t",
  row.names = F,
  col.names = T
)

################################################################################
# Plot WLM results

## Maternal DUO estimates

duos <- c("Maternal (fetal adjusted)", "Maternal (unadjusted)")

mr_res <- fread("wlm_500_0.2_mr_analysis.txt") %>%
  left_join(., out_info, by = c("outcome" = "outcome")) %>% ## annotate outcomes
  filter(adjustment %in% duos) %>%
  mutate(
    exposure = factor(gene.exposure, levels = exp_info$gene_name),
    adjustment = factor(adjustment,
    levels = duos)
  )

colours <- c('#8C103D', '#D81B60', '#FE6100', '#FFC107')
names(colours) <- unique(mr_res$adjustment)

plot_mr_special(mr_res, "adjustment") %>%
  ggsave(
    filename = "lpa_wlm_mr_res_mat_duos.png",
    width = 9,
    height = 7,
    dpi = 300)

# TRIO estimates

trios <- c(
  "Paternal (maternal+fetal adjusted)", "Paternal (unadjusted)",
  "Fetal (maternal+paternal adjusted)", "Fetal (maternal adjusted)",
  "Fetal (unadjusted)", "Maternal (fetal+paternal adjusted)",
  "Maternal (fetal adjusted)", "Maternal (unadjusted)"
  )

mr_res <- fread("wlm_500_0.2_mr_analysis.txt") %>%
## annotate outcomes
  left_join(., out_info, by = c("outcome" = "outcome")) %>% 
  mutate(
    exposure = factor(gene.exposure, levels = exp_info$gene_name),
    adjustment = factor(adjustment,
    levels = trios
  )

colours <- c(
  '#8C103D', '#D81B60', '#FF5B97', '#FE6100',
  '#FFC107', '#FFE749', '#125D9E', '#1E88E5'
  )

names(colours) <- unique(mr_res$adjustment)

plot_mr_special(mr_res, "adjustment") %>%
  ggsave(
    filename = "lpa_wlm_mr_res_mat_trios.png",
    width = 9,
    height = 9,
    dpi = 300)
