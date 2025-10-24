################################################################################
#                        Genetic score LPA x APPOs in UK Biobank                         #
################################################################################

##
# This script creates a Lp(a) genetic score (GS) for the MR-PREG UKB sample
# based on weights from the GWAS by Mack et al (2017) and tests its association 
# with APPOs across quintiles of the GS distribution, as well as overall.
##

# N.B. To run this script, run:
# system("dx download scripts/UKB_RAP/3_gs_lpa_appos.R")
# and then
# source("3_gs_lpa_appos.R")
# in the RStudio console of your UKB-RAP Posit session.

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

## Download required data to RStudio session

# Phenotypic dataset
system('dx download -f data/UKB_RAP/lpa_appos_dataset.tsv', intern = TRUE)

# GS SNP list and weights
system('dx download -f data/UKB_RAP/lpa_mack2017_weights.txt', intern = TRUE)

# APPO info
system('dx download -f scripts/info.xlsx', intern = TRUE)

# Helper functions
system('dx download -f scripts/UKB_RAP/RAP_functions.R', intern = TRUE)

################################################################################
# Build analysis dataset

source("RAP_functions.R")

## Read in data
out_info <- readxl::read_excel("info.xlsx")
dataset <- readr::read_tsv("lpa_appos_dataset.tsv")
varlist_gs <- readr::read_tsv("lpa_mack2017_weights.txt")

## Format variant data for GS function
vars <- varlist_gs %>%
  dplyr::select(
	rsid = SNP,
	effect_allele,
	other_allele,
	beta,
  chr,
  pos
	) 

## Create GS using ukbrapR
lpa_gs <- ukbrapR::create_pgs(
  in_file = vars,
  pgs_name = "lpa_gs",
  out_file = "lpa_pgs_output",
  source = "imputed",
  progress = TRUE,
  very_verbose = TRUE,
  # overwrite = TRUE
  )

## Merge GS with dataset
gs_dataset <- dataset %>%
  dplyr::left_join(., lpa_gs, by = "eid") %>%
# Convert GS to SD units
  dplyr::mutate(
    gs = scale(lpa_gs, )[,1],
# Generate GS quintiles 
    gs_quint = ntile(gs, 5)
  )

## Visualise distribution of GS
png("lpa_gs_distribution.png", width = 800, height = 600)
hist(gs_dataset$gs, breaks = 100)
dev.off()
summary(gs_dataset$gs)

## Visualise GS quintiles
table(gs_dataset$gs_quint)

dim(gs_dataset)
# [1] 224313     58

## Save dataset with GS
readr::write_tsv(gs_dataset, "lpa_appos_dataset_with_pgs.tsv")
system('dx upload lpa_appos_dataset_with_pgs.tsv --path data/UKB_RAP/', intern = TRUE)

## Save GS distribution plot
system('dx upload lpa_gs_distribution.png --path results/genetic_score/', intern = TRUE)

################################################################################
# Look at event numbers by GS quintile

## Create counts of events in out_info and venous events per quintile

outcome_list <- c(
  "any_event", 
  "venous_event", 
  "atheros_event", 
  out_info$outcome[out_info$cat == "bin"]
  )

count_list <- list()

## Loop through each outcome
for (outcome in outcome_list) {

  if (outcome %in% colnames(gs_dataset)) {
    subset_df <- gs_dataset[gs_dataset[[outcome]] == 1, ]
    # Count cases per quintile
    quintile_counts <- table(factor(subset_df$gs_quint, levels = 1:5))
  } else {
    subset_df <- data.frame()
    quintile_counts <- rep(NA, 5)
  }
  
  # Create a df row
  count_row <- data.frame(
    outcome = outcome,
    quint_1 = quintile_counts[1],
    quint_2 = quintile_counts[2],
    quint_3 = quintile_counts[3],
    quint_4 = quintile_counts[4],
    quint_5 = quintile_counts[5]
  )
  
  # Append to count list
  count_list[[outcome]] <- count_row
}

## Combine all rows into a single data frame
count_df <- do.call(rbind, count_list) %>%
  filter(!is.na(quint_1))  # Remove rows with NA counts (outcomes not in dataset)

## View the counts
print(count_df)

## Save counts to a file
readr::write_tsv(count_df, "event_counts_by_gs_quintile.tsv")
system('dx upload event_counts_by_gs_quintile.tsv --path results/genetic_score/', intern = TRUE)

################################################################################
# Looking at Lp(a) GS x Lp(a) levels

## Test association between Lp(a) genetic score and Lp(a) levels
pcs <- paste0("PC_", 1:40)
geno_batch <- "geno_batch"
lpa_test <- lm(reformulate(c("gs", pcs, geno_batch), response = "lpa_mes"), data = gs_dataset)
coef(summary(lpa_test))

## Save results
lpa_test_df <- broom::tidy(lpa_test, conf.int = TRUE) %>%
  select(estimate, std.error, p.value, conf.low, conf.high)

readr::write_tsv(lpa_test_df, "lpa_gs_x_lpa.tsv")
system('dx upload lpa_gs_x_lpa.tsv --path results/genetic_score/', intern = TRUE)

## Plot violin plot of lp(a) levels across GS quintiles
colours <- c(
  "1" = '#0016deff',
  "2" = '#77abffff',
  "3" = '#bd9eeaff',
  "4" = '#ffa6beff',
  "5" = '#ff5454ff'
)

medians <- gs_dataset %>%
    group_by(gs_quint) %>%
    summarise(
      median_lpa = median(lpa_mes, na.rm = TRUE),
      .groups = "drop"
    )

png("lpa_levels_by_gs_quintile.png", width = 12, height = 8, units = "in", res = 300)
ggplot(gs_dataset, aes(x = factor(gs_quint), y = lpa_mes, fill = factor(gs_quint))) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Add violin plot with 50% transparency
  geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.7, outlier.shape = NA) +  # Hide outlier points
  scale_fill_manual(values = colours) +  # Set colors for the violins
  scale_y_continuous(limits = c(0, NA)) +  # Set lower y-axis limit to 0, keep upper limit automatic
  # Add median values as text labels on top
  geom_text(
    data = medians,
    aes(x = gs_quint, y = median_lpa, label = sprintf("%.1f", median_lpa)),
    vjust = -0.5,
    size = 3.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Lp(a) Genetic Score Quintile",
    y = "Lp(a) Levels (nmol/L)",
    caption = "N per quintile = 35097"
  ) +
  theme_minimal() +
  guides(fill = "none")
dev.off()

system('dx upload lpa_levels_by_gs_quintile.png --path results/genetic_score/', intern = TRUE)

## Plot histogram of Lp(a) levels by GS quintile - 5 separate rows 

png("lpa_levels_hist_by_gs_quintile.png", width = 16, height = 12, units = "in", res = 300)
ggplot(gs_dataset, aes(x = lpa_mes, fill = factor(gs_quint))) +
  geom_histogram(bins = 50, alpha = 0.8, color = "white", size = 0.2) +
  facet_wrap(~ paste("Quintile", gs_quint), ncol = 1, scales = "free_y", 
             labeller = labeller(.cols = function(x) paste("Genetic Score", x))) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(
    x = "Lp(a) Levels (nmol/L)",
    y = "Count",
    fill = "GS Quintile"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.position = "none"
  )
dev.off()

system('dx upload lpa_levels_hist_by_gs_quintile.png --path results/genetic_score/', intern = TRUE)

################################################################################
# Association testing - overall and quintile-stratified analyses

## Create dataframe for function arguments to iterate through
args <- select(out_info, outcome, cat)

## Add the composite outcomes to the args dataframe
composite_outcomes <- data.frame(
  outcome = c("any_event", "venous_event", "atheros_event"),
  cat = rep("bin", 3)
)
args <- bind_rows(args, composite_outcomes) %>%
  # Only keep outcomes that exist in dataset
  filter(outcome %in% colnames(gs_dataset))  

## Overall analysis - for comparison
overall_results <- purrr::pmap_dfr(
  args,
  assoc_test,
  df = gs_dataset,
  score_var = "gs",
  quint_var = NULL 
)

## Quintile-stratified analysis
quintile_results <- purrr::pmap_dfr(
  args,
  assoc_test,
  df = gs_dataset,
  score_var = "gs",
  quint_var = "gs_quint"
)

## Combine results
combined <- rbind(
  overall_results %>% mutate(quintile = "overall"),
  quintile_results
) %>% mutate(FDR = p.adjust(p_value, method = "fdr"))

## Save results
readr::write_tsv(combined, "lpa_gs_x_appo.tsv")
system('dx upload lpa_gs_x_appo.tsv --path results/genetic_score/', intern = TRUE)

################################################################################
# Plot results

to_plot <- combined %>%
  left_join(., out_info, by = "outcome") %>%
  mutate(
    outcome_full = case_when(
      outcome == "any_event" ~ "Any thrombotic event",
      outcome == "venous_event" ~ "Venous thrombotic event",
      #outcome == "atheros_event" ~ "Atherosclerotic thrombotic event",
      TRUE ~ out_info$outcome_full[match(outcome, out_info$outcome)]
    ),
    cat = case_when(
      outcome %in% c("any_event", "venous_event", "atheros_event") ~ "bin",
      TRUE ~ out_info$cat[match(outcome, out_info$outcome)]
    ),
    quintile = as.character(quintile),
    b = -b # reverse direction for comparability with rest of study
  ) %>%
  filter(!is.na(outcome_full))

## Define colours - required for plot function to run
colours <- c(
  "overall" = '#000000ff',
  "1" = '#0016deff',
  "2" = '#77abffff',
  "3" = '#bd9eeaff',
  "4" = '#ffa6beff',
  "5" = '#ff5454ff'
)

## Plot 
png("lpa_gs_x_appo.png", width = 14, height = 13, units = "in", res = 300)
plot_strat(to_plot, "quintile")
dev.off()

## Plot without thrombotic events
png("lpa_gs_x_appo_no_thrombotic.png", width = 14, height = 12, units = "in", res = 300)
plot_strat(
  to_plot %>% filter(!outcome %in% c("any_event", "venous_event", "atheros_event")), 
  "quintile"
  )
dev.off()

## Upload plots
system('dx upload lpa_gs_x_appo.png --path results/genetic_score/', intern = TRUE)
system('dx upload lpa_gs_x_appo_no_thrombotic.png --path results/genetic_score/', intern = TRUE)