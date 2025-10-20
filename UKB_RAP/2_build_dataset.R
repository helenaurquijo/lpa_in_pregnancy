################################################################################
#                        PRS LPA x APPOs in UK Biobank                         #
################################################################################

##
# This script builds a dataset from 1) the derived MR-PREG phenotypes and
# 2) ICD9/10 main diagnoses for thrombotic events, 3) lipoprotein(a) levels and
# 4) 40 genetic principal components, for further analysis.
##

# to upload this script to DNAnexus, run:
# system("dx upload scripts/UKB_RAP/2_build_dataset.R")

## Load packages
if(!require(pacman)) install.packages("pacman")
pacman::p_load(dplyr, readr, readxl, stringr)

## Download required data to RStudio session

# MR-PREG phenotypes
system('dx download -f data/derived/Phenotype_MRPREG.txt', intern = TRUE)
system('dx download -f data/derived/phenotypic_data_maternity.txt', intern = TRUE)

# eid x IEU ID linker file
system('dx download -f data/derived/linker_app2393.csv', intern = TRUE)

# ICD9/10 main diagnoses and lp(a)
system('dx download -f data/extracted_data.tsv', intern = TRUE)

# APPO info
system('dx download -f scripts/info.xlsx', intern = TRUE)

# Helper functions
system('dx download -f scripts/RAP_functions.R', intern = TRUE)

## Read in data
derived_data <- readr::read_tsv("Phenotype_MRPREG.txt")
derived_data_mat <- readr::read_tsv("phenotypic_data_maternity.txt")
linker <- readr::read_csv("linker_app2393.csv") %>% rename(eid = app)
out_info <- readxl::read_excel("info.xlsx")
extdat <- readr::read_tsv("extracted_data.tsv")
source("RAP_functions.R")

## Merge maternity and questionnaire-based derived data

dim(derived_data)
# [1] 224313     73
dim(derived_data_mat)

## join keeping all rows from both datasets
all_derived <- full_join(
  derived_data,
  derived_data_mat,
  by = "FID"
)

dim(all_derived)

## Match MR-PREG data with participant IDs using linker

dataset <- all_derived %>%
  dplyr::select(FID, any_of(out_info$outcome), smoking, bmi_kgm2) %>%
  dplyr::left_join(., linker, by = c("FID" = "ieu"))

dim(dataset)
# [1] 224313     12

## Merge derived dataset with additional data

dataset1 <- dataset %>%
  dplyr::left_join(., extdat, by = "eid") %>%
  rename_with(
    # Rename p22009_a1-40 to PC_1-40 (only if columns exist)
    ~ paste0("PC_", 1:40),
    any_of(paste0("p22009_a", 1:40))
  ) %>%
  rename(
    ICD10 = 'p41202',      # HES main ICD-10 codes
    ICD9 = 'p41203',       # HES main ICD9 codes
    #BMI = 'p21001_i0',      # Body mass index
    #pack_years_smoking = 	'p20162_i0', # Pack years adult smoking as proportion of life span exposed to smoking
    lpa_mes = 'p30790_i0', # Lipoprotein(a) measurements
    exclude = 'p22010',  # Recommended genetic analysis exclusions
    geno_batch = 'p22000', # Genotyping batch
    geno_ethnic = 'p22006',  # Genetic ethnic grouping
    related_excl = 'p22018',  # Genetic relatedness exclusions
    aneuploidy = 'p22019',  # Sex chromosome aneuploidy
    used_in_pc = 'p22020',  # Used in genetic principal components - a (maximal) subset of unrelated participants after applying some QC filtering as described in the Supplement of Bycroft et al., 2018.
  )

## Peform initial data QC

dataset2 <- dataset1 %>%
  dplyr::filter(
    is.na(exclude) &
    is.na(related_excl) &
    is.na(aneuploidy) &
    geno_ethnic == 1 &
    used_in_pc == 1
  )

## Prune related individuals

## Get thrombotic events

dataset3 <- dataset2 %>%
  dplyr::mutate(
    # Handle missing ICD columns
    ICD10 = tidyr::replace_na(ICD10, ""),
    ICD9 = tidyr::replace_na(ICD9, ""),
    
    # ICD10 venous thrombotic event codes (ICD10)
    icd10_venous = str_detect(ICD10, paste(c(
      "O222",  # Superficial thrombophlebitis in pregnancy
      "O223",  # Deep phlebothrombosis in pregnancy
      "O225",  # Cerebral venous thrombosis in pregnancy
      "O229",  # Venous complication in pregnancy, unspecified
      "O870",  # Superficial thrombophlebitis in the puerperium
      "O871",  # Deep phlebothrombosis in the puerperium
      "O873",  # Cerebral venous thrombosis in the puerperium
      "O879",  # Venous complication in the puerperium, unspecified
      "O882"   # Obstetric thromboembolism
    ), collapse = "|")),
    
    # ICD10 atherosclerotic event codes (ICD10)
    icd10_atheros = str_detect(ICD10, "0994"),
    
    # ICD9 venous thrombotic event codes (ICD9)
    icd9_venous = str_detect(ICD9, paste(c(
      "6712",  # Superficial thrombophlebitis in pregnancy and puerperium
      "6713",  # Deep phlebothrombosis antepartum
      "6715",  # Other phlebitis and thrombosis
      "6719",  # Unspecified venous complication
      "6714",  # Deep phlebothrombosis postpartum
      "6732",  # Obstetrical blood-clot embolism
      "6740"   # Cerebrovascular disorders in the puerperium
    ), collapse = "|")),
    
    # ICD9 atherosclerotic event codes (ICD9)
    icd9_atheros = str_detect(ICD9, "6486"),

    # Combined event indicators (handle NAs)
    venous_event = ifelse(icd10_venous | icd9_venous, 1, 0),
    atheros_event = ifelse(icd10_atheros | icd9_atheros, 1, 0),
    any_event = ifelse(venous_event == 1 | atheros_event == 1, 1, 0)
  )

dim(dataset3)
# [1] 224313     67


## Create unified tally of all variables and events 

tally_data <- list()

# 1. ICD10 venous codes
icd10_venous_codes <- c("O222", "O223", "O225", "O229", "O870", "O871", "O873", "O879", "O882")
icd10_descriptions <- c(
  "Superficial thrombophlebitis in pregnancy", "Deep phlebothrombosis in pregnancy",
  "Cerebral venous thrombosis in pregnancy", "Venous complication in pregnancy, unspecified",
  "Superficial thrombophlebitis in the puerperium", "Deep phlebothrombosis in the puerperium",
  "Cerebral venous thrombosis in the puerperium", "Venous complication in the puerperium, unspecified",
  "Obstetric thromboembolism"
)
tally_data[["icd10"]] <- data.frame(
  Category = rep("ICD10 Venous", length(icd10_venous_codes)),
  Code = icd10_venous_codes,
  Description = icd10_descriptions,
  Count = sapply(icd10_venous_codes, function(code) {
    sum(str_detect(dataset3$ICD10, code), na.rm = TRUE)
  }),
  stringsAsFactors = FALSE
)

# 2. ICD9 venous codes
icd9_venous_codes <- c("6712", "6713", "6715", "6719", "6714", "6732", "6740")
icd9_descriptions <- c(
  "Superficial thrombophlebitis in pregnancy and puerperium", "Deep phlebothrombosis antepartum",
  "Other phlebitis and thrombosis", "Unspecified venous complication",
  "Deep phlebothrombosis postpartum", "Obstetrical blood-clot embolism",
  "Cerebrovascular disorders in the puerperium"
)
tally_data[["icd9"]] <- data.frame(
  Category = rep("ICD9 Venous", length(icd9_venous_codes)),
  Code = icd9_venous_codes,
  Description = icd9_descriptions,
  Count = sapply(icd9_venous_codes, function(code) {
    sum(str_detect(dataset3$ICD9, code), na.rm = TRUE)
  }),
  stringsAsFactors = FALSE
)

# 3. Event indicators
event_indicators <- c("any_event", "venous_event", "atheros_event")
event_descriptions <- c(
  "Any thrombotic event (venous or atherosclerotic)",
  "Venous thrombotic event",
  "Atherosclerotic event"
)
tally_data[["events"]] <- data.frame(
  Category = rep("Event Indicators", length(event_indicators)),
  Code = event_indicators,
  Description = event_descriptions,
  Count = sapply(event_indicators, function(ind) {
    count <- sum(dataset3[[ind]], na.rm = TRUE)
    sprintf("%d (%.2f%%)", count, count/nrow(dataset3)*100)
  }),
  stringsAsFactors = FALSE
)

# 4. Outcome variables
if(length(out_info$outcome) > 0) {
  outcomes_present <- intersect(out_info$outcome, names(dataset3))
  if(length(outcomes_present) > 0) {
    # Create mapping of outcomes to types and descriptions
    outcome_types <- setNames(out_info$cat, out_info$outcome)
    outcome_descriptions <- setNames(out_info$outcome_full, out_info$outcome)
    
    tally_data[["outcomes"]] <- data.frame(
      Category = rep("Outcomes", length(outcomes_present)),
      Code = outcomes_present,
      Description = sapply(outcomes_present, function(o) {
        ifelse(!is.na(outcome_descriptions[o]), outcome_descriptions[o], o)
      }),
      Count = sapply(outcomes_present, function(o) {
        summarize_variable(dataset3, o, outcome_types[o])
      }),
      stringsAsFactors = FALSE
    )
  }
}

# Combine all tallies
all_tallies <- do.call(rbind, tally_data)
rownames(all_tallies) <- NULL

# Display unified tally
print("Unified tally of all variables:")
print(all_tallies)

# Save the tally to a file
write.csv(all_tallies, "variable_tally.csv", row.names = FALSE)
system('dx upload variable_tally.csv --dest data/', intern = TRUE)

## Clean analysis dataset

final_dataset <- dataset3 %>%
  dplyr::select(
    eid,
    any_of(out_info$outcome),  
    any_event,
    venous_event,
    atheros_event,
    bmi_kgm2,
    smoking,
    lpa_mes,
    geno_batch,
    any_of(paste0("PC_", 1:40))  
  )

dim(final_dataset)
#[1] 224313     55

## Save dataset
readr::write_tsv(final_dataset, "lpa_appos_dataset.tsv")
system('dx upload lpa_appos_dataset.tsv --dest data/', intern = TRUE)