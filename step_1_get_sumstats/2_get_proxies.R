################################################################################
#                          Get proxies for Lp(a) IVs                           #
################################################################################

##
# This script aims to extract a list of tagging variants for each of the 
# SNPs being used as IVs as a lot of them are missing from the outcome datasets
#
# Went for this approach as opposed to using plink 1.9 as the --show-tags
# function does not provide the tag allele in phase with index allele
#
# This script is based off of code by Lizzy Aiton which is on Github: 
# eaiton/autoimmune-pregnancy/blob/main/autoimmune/2_a_sel_proxies.R
##

print(Sys.time())

library(data.table)
library(dplyr)
library(readxl)
library(TwoSampleMR)
library(MRutils)

## Load env and paths
readRenviron("scripts/config.env")
projectdir <- Sys.getenv("projectdir")

proxy_cache <- paste0(projectdir, "/data/exposures/proxy_cache_EUR/")

###############################################################################
# Find proxies for all exposure instruments

setwd(paste0(projectdir, "/data/exposures/instruments"))

instr <- fread("rsids_main_lpa.txt", header = F)

proxies <- MRutils::get_proxies(
  rsids = instr$V1,
  token = Sys.getenv("LDLINK_TOKEN"),
    # request a personal LDlink token here:
    # https://ldlink.nih.gov/?tab=apiaccess
  population = "EUR",
  results_dir = proxy_cache,
  skip_api = FALSE,
  r2_threshold = 0.8)
    
head(proxies)
dim(proxies)

## Save proxy info
fname <- "proxies_info_main_lpa.txt"
fwrite(proxies, fname, sep ="\t", col.names = TRUE)

## Save list of proxy SNPs
fname <- "proxies_rsid_main_lpa.txt"
rsids <- as.data.frame(unique(proxies$rsid))
fwrite(rsids, fname, sep ="\t", col.names = FALSE)

print("Potential proxy SNPs saved")

q("no")
