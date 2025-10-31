################################################################################
#                        LoF LPA x APPOs in UK Biobank                         #
################################################################################

##
# 
##

# N.B. To run this script, run:
# system("dx download scripts/UKB_RAP/4_lof_lpa_appos.R")
# and then
# source("4_lof_lpa_appos.R")
# in the RStudio console of your UKB-RAP Posit session. Say yes to creating a
# personal library and updating all packages when prompted.

################################################################################
# Set up environment

## Load packages
# Ensure pacman is installed
if (!require("pacman")) install.packages("pacman")

# Force install latest ggplot2 from CRAN to avoid compatibility issues
install.packages("ggplot2", dependencies = NA)

# Load rest of packages using pacman
pacman::p_load(
  rlang, dplyr, readr, readxl, stringr, purrr, broom,
  ggplot2, cowplot
)

## Install ukbrapR and ggforesplot packages from GitHub
remotes::install_github("lcpilling/ukbrapR")
library(ukbrapR)

devtools::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)
#  N.B ggforest uses a deprecated ggplot2 function, using suppressWarnings()
#     until this is fixed. 

## Install Bioconductor utilities
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

## Load data onto RStudio session

# Phenotypic dataset
system('dx download -f data/UKB_RAP/lpa_appos_dataset.tsv', intern = TRUE)
dataset <- readr::read_tsv("lpa_appos_dataset.tsv")

# APPO info
system('dx download -f scripts/info.xlsx', intern = TRUE)
out_info <- readxl::read_excel("info.xlsx")

# Helper functions
system('dx download -f scripts/UKB_RAP/RAP_functions.R', intern = TRUE)
source("RAP_functions.R")


################################################################################
# Prepare LoF variant list

## List of LoF variants in LPA
#   These were found in 10.1016/j.atherosclerosis.2022.04.003 and GeneBass
vars <- c(
   "rs41259144",
   "rs41272114",
   "rs41267813",
   "rs41267811",
   "rs143431368",
   "rs762351463",
   "rs41267811",
   "rs200099994",
   "rs559580002"
)

## Get locations of SNPs
mart <- useEnsembl(
  biomart = "snps",
  )

mart <- useMart(biomart="ENSEMBL_MART_SNP", host="https://grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_snp")

lof_list <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start"),
  filters = "snp_filter",
  values = vars,
  mart = mart
)

varlist <- lof_list %>%
  rename(
    chr = chr_name,
    pos = chrom_start,
    rsid = refsnp_id
  )

################################################################################
# Extract genotypes from UKB RAP

## Extract genotypes
genotypes <- extract_variants(
  varlist,
  out_bed = "lpa_lof",
  source = "imputed",
  progress = TRUE,
  very_verbose = TRUE
  )

## Merge genotypes with dataset
lof_dataset <- dataset %>%
  inner_join(genotypes, by = "eid")

print(dim(lof_dataset))
# [1] 175485     71

## Create lof_homozygous variable
lof_dataset <- lof_dataset %>%
  mutate(
    lof_hom = if_else(
      rowSums(across(all_of(colnames(genotypes)[-1]), ~ . == 2), na.rm = TRUE) > 0,
      1,
      0
    )
  )

table(lof_dataset$lof_hom)
#      0      1 
# 175168    317 

## Save updated dataset
write_tsv(lof_dataset, "lpa_appos_dataset_with_lof.tsv")
system('dx upload lpa_appos_dataset_with_lof.tsv --path data/UKB_RAP/', intern = TRUE)

################################################################################
# Descriptive stats and plots

## Create table with LoF variant counts
variant_counts <- data.frame(
  variant = character(),
  WT = integer(),
  HET = integer(),
  HOM = integer(),
  stringsAsFactors = FALSE
)

for (v in colnames(genotypes)[-1]) {
  counts <- table(lof_dataset[[v]], useNA = "ifany")
  
  # Convert table to named vector for easier access
  counts_vector <- as.vector(counts)
  names(counts_vector) <- names(counts)
  
  # Extract counts (handling potential missing genotypes)
  wt_count <- counts_vector["0"]
  het_count <- counts_vector["1"]
  hom_count <- counts_vector["2"]
  
  # Handle NAs by setting to 0 if category doesn't exist
  if (is.na(wt_count)) wt_count <- 0
  if (is.na(het_count)) het_count <- 0
  if (is.na(hom_count)) hom_count <- 0
  
  # Add to table
  variant_counts <- rbind(variant_counts, data.frame(
    variant = v,
    WT = wt_count,
    HET = het_count,
    HOM = hom_count
  ))
}

## Print the summary table
print(variant_counts)

## Save the table
write.csv(variant_counts, file = "lof_variant_counts.csv", row.names = FALSE)
system('dx upload lof_variant_counts.csv --path results/lof/', intern = TRUE)

## Compare case numbers: lof homozygotes vs everyone else
outcome_list <- c(
  "any_event", 
  "venous_event", 
  "atheros_event",
  "cad_event",
  out_info$outcome[out_info$cat == "bin"]
  )

count_list <- list()

## Loop through each outcome
for (outcome in outcome_list) {

  if (outcome %in% colnames(lof_dataset)) {
    subset_df <- lof_dataset[lof_dataset[[outcome]] == 1, ]
    # Count cases per quintile
    counts <- table(factor(subset_df$lof_hom, levels = c(0, 1)))
  } else {
    subset_df <- data.frame()
    quintile_counts <- rep(NA, 5)
  }

  # Create a df row
  count_row <- data.frame(
    outcome = outcome,
    lof_hom = counts[1],
    lof_het_wt = counts[2]
  )
  
  # Append to count list
  count_list[[outcome]] <- count_row
}

## Combine all rows into a single data frame
count_df <- do.call(rbind, count_list)

## View the counts
print(count_df)

## Save counts to a file
readr::write_tsv(count_df, "event_counts_by_lof_geno.tsv")
system('dx upload event_counts_by_lof_geno.tsv --path results/lof/', intern = TRUE)

## Look at lpa levels by variant - violin plot per variant - comparing WT/HET/HOM

## Create individual plots for each variant
variant_plots <- list()
for (v in colnames(genotypes)[-1]) {
  # Skip variants with no variation
  if (length(unique(na.omit(lof_dataset[[v]]))) <= 1) {
    cat("Skipping", v, "- insufficient variation\n")
    next
  }
  
  # Create plot
  plot <- create_violin_plot(lof_dataset, v)
  variant_plots[[v]] <- plot
  
  # Save individual plot
  ggsave(
    filename = paste0(v, "_lpa_violin.png"),
    plot = plot,
    width = 7,
    height = 5,
    dpi = 300
  )
}

## Download all plots
for (v in names(variant_plots)) {
  system(paste0('dx upload ', v, '_lpa_violin.png --path results/lof/plots/'), intern = TRUE)
}

## Compare number of individuals out of reportable range by genotype
reportable_counts <- dplyr::select(lof_dataset, lof_hom, lpa_mes, lpa_report) %>%
  mutate(
    reportable = case_when(
      lpa_report == 4 ~ "too_low",
      lpa_report == 5 ~ "too_high",
      lpa_report == 1 ~ "within_range",
      TRUE ~ "NA"
    )
  ) %>%
  group_by(lof_hom, reportable) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = reportable, values_from = n, values_fill = 0)

print(reportable_counts)

# Save reportable counts
readr::write_tsv(reportable_counts, "lpa_report_by_lof_geno.tsv")
system('dx upload lpa_report_by_lof_geno.tsv --path results/lof/', intern = TRUE)

################################################################################
# Test association of homozygosity for LoF variants with APPOs

## Create dataframe for function arguments to iterate through
args <- dplyr::select(out_info, outcome, cat)

## Add the composite outcomes to the args dataframe
composite_outcomes <- data.frame(
  outcome = c("any_event", "venous_event", "atheros_event", "cad_event"),
  cat = rep("bin", 4)
)
args <- bind_rows(args, composite_outcomes) %>%
  # Only keep outcomes that exist in dataset
  filter(outcome %in% colnames(lof_dataset))

## Run association tests - reutilising function from GS analysis
results <- purrr::pmap_dfr(
  args,
  assoc_test,
  df = lof_dataset,
  score_var = "lof_hom",
  quint_var = NULL # No quintiles for this analysis
) %>% 
  dplyr::select(-quintile) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))

## Save results
readr::write_tsv(results, "lof_hom_x_appos.tsv")
system('dx upload lof_hom_x_appos.tsv --path results/lof/', intern = TRUE)

###################################################################################
# Plot results

## Format results
to_plot <- results %>%
  left_join(., out_info, by = "outcome") %>%
  mutate(
    outcome_full = case_when(
      outcome == "any_event" ~ "Any thrombotic event",
      outcome == "venous_event" ~ "Venous thrombotic event",
      #outcome == "atheros_event" ~ "Atherosclerotic thrombotic event",
      outcome == "cad_event" ~ "Coronary artery disease",
      TRUE ~ out_info$outcome_full[match(outcome, out_info$outcome)]
    ),
    cat = case_when(
      !outcome %in% out_info$outcome ~ "bin",
      TRUE ~ out_info$cat[match(outcome, out_info$outcome)]
    )
  ) %>%
  filter(!is.na(outcome_full)) %>%
  # remove extreme estimates for plotting purposes
  filter(abs(b) < 10)


## Plot 
suppressWarnings(
  cont <- to_plot %>%
    filter(cat == "cont") %>%
    forestplot(
      df = .,
      name = outcome_full,
      estimate = b,
      logodds = FALSE,
      xlab = "Mean difference in trait (95% CI)",
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      legend.position = "none"
    ) + ggtitle("Continuous outcomes") +
    scale_color_manual(values = colours)
)

suppressWarnings(
  bin <- to_plot %>%
    filter(cat == "bin") %>%
    forestplot(
      df = .,
      name = outcome_full,
      estimate = b,
      logodds = TRUE,
      xlab = "Odds ratio for APPO (95% CI)",
      pvalue = FDR
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank()
    ) + ggtitle("Binary outcomes") +
    scale_color_manual(values = colours)
)

png("lpa_lof_x_appo.png", width = 14, height = 12, units = "in", res = 300)  

suppressWarnings(
  cowplot::plot_grid(
    bin, cont, 
    ncol = 1, 
    rel_heights = c(6, 1),
    align = "v",
    labels = NULL
  ) %>%
    cowplot::plot_grid(., ncol = 1, rel_heights = c(0.1, 1))
)

dev.off()

system('dx upload lpa_lof_x_appo.png --path results/lof/', intern = TRUE)
