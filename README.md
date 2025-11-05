# lpa_in_pregnancy
Supporting information and scripts for the study of the effects of the modulation of lp(a) levels on adverse pregnancy and perinatal outcomes (APPOs). UoB x NVS collaboration.

## Workspace Setup

### 1. Create Directory Structure
Run the setup script to create the required folder structure:
```bash
bash step_0_set_up.sh
```
which should look like:

```
.
|-- data
|   |-- UKB_RAP
|   |-- exposures
|   |   |-- full_sumstats
|   |   |-- instruments
|   |   `-- raw
|   |-- for_mrlap
|   `-- outcomes
|       |-- full_sumstats
|       `-- raw
|-- logs
|-- results
|   |-- UKB_RAP
|   |   |-- genetic_score
|   |   `-- lof
|   |       `-- plots
|   |-- iv_checks
|   |-- main_analysis
|   `-- sensitivity_analyses
|       `-- assoc_plots
`-- scripts
    |-- UKB_RAP
    |-- functions
    |-- step_1_get_sumstats
    |-- step_2_run_mr
    |-- step_3_other_sens
    `-- step_4_iv_checks
```

with the **contents of this repository located in the `scripts` folder**.

### 2. Set up aliases
Create a `config.env` file in the scripts directory using the template `config_template.env`. N.B. you must set up access tokens for [OpenGWAS](https://opengwas.io) and [LDLink](https://ldlink.nih.gov/?tab=apiaccess).

### 3. Software dependencies
This project was run in R (v4.3.3) using an [renv](https://rstudio.github.io/renv/) environment. To install the required packages, open R in the scripts folder and:
1. Install renv using 
```R
install.packages("renv")
```
2. Install the relevant packages using 
```R
source("renv/activate.R")
renv::restore() # this step will take a while
```
### 4. Data required
Besides access to the MR-PREG meta-analysis and study-specific GWAS full summary statistics, the following datasets have to be downloaded:

- Reference panel of 10K UKB participants (private dataset but easily replaced with 1KG3)
- Inputs for LDSC as part of MRlap, accessed [here](https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v).
- Olink [UKB-PPP](https://www.nature.com/articles/s41586-023-06592-6) full summary statistics for LPA and metadata files, accessed [here](https://www.synapse.org/Synapse:syn51364943). 
- Lead lp(a) GWAS summary statistics from [Mack et al., 2017](https://linkinghub.elsevier.com/retrieve/pii/S0022227520339547), accessed [here](https://www.jlr.org/cms/10.1194/jlr.M076232/attachment/cfbb7466-d752-4ab7-b82b-558b638259e2/mmc1.pdf). These had to be extracted manually from a PDF file but the steps are documented in `1_get_ivs.R`.
- PPH full GWAS summary statistics from [Westergaard et al., 2024](https://www.nature.com/articles/s41588-024-01839-y), accessed [here](https://download.decode.com/form/2024-steinthorsdottir-postpartu).
- Lp(a) full GWAS summary statistics from [Sinnott-Armstrong et al., 2021](https://www.nature.com/articles/s41588-020-00757-z), accessed [here](https://nih.figshare.com/articles/dataset/The_meta-analyzed_GWAS_summary_statistics_for_35_lab_biomarkers_described_in_Genetics_of_35_blood_and_urine_biomarkers_in_the_UK_Biobank_/12355382?file=22768706).

Any other data is downloaded by this repository's scripts.  

## Analysis Workflow
The analysis is organized into sequential steps:

### Step 1: Get Summary Statistics (`step_1_get_sumstats/`)
- `1_get_ivs.R`: Extract instrumental variables (IVs)
- `2_get_proxies.R`: Find proxy SNPs for IVs
- `3_out_SNP_extract.sh`: Extract outcome SNPs from MR-PREG
- `4_out_SNP_extract.R`: Process extracted SNPs
- `5_out_SNP_extract_other.R`: Extract SNPs for non-MR-PREG outcomes
- `6_get_full_sumstats.R`: Extract summary statistics from Lp(a), PPH and conmal_any GWAS
- `7_out_full_extract.sh`: Extract full MR-PREG summary statistics

### Step 2: Run MR Analyses (gIVW method) (`step_2_run_mr/`)
- `8_main_MR.R`: Main Mendelian randomization analysis
- `9_alt_lpa_MR.R`: Alternative Lp(a) exposure datasets
- `10_wlm_MR.R`: Weighted linear model - adjusting for fetal genetic effects
- `11_loo_MR.R`: Leave-one-out of MR-PREG studies

### Step 3: Other Sensitivity Analyses (`step_3_other_sens/`)
- `12_run_mrlap.R`: Run MRlap for sample overlap correction
- `13_coloc.R`: Genetic colocalisation: Lp(a) levels x APPOs

### Step 4: IV Validation (`step_4_iv_checks/`)
- `14_pos_ctrl_MR.R`: Positive control analysis

### UKB-RAP analysis (`UKB_RAP/`)
The scripts in this folder were run in the [UKB RAP](https://www.ukbiobank.ac.uk/use-our-data/research-analysis-platform/).

- `1_data_extraction.ipynb`: Extraction of relevant data fields from the database. Run in a Spark Cluster instance of Jupyter lab.
- `2_build_dataset.R`: Merging of data with MR-PREG dataset and formatting. Run on an RStudio Workbench session.
- `3_gs_lpa_appos.R`: Investigating the association between Lp(a) genetic score and APPOs. Run on an RStudio Workbench session.
- `4_lof_lpa_appos.R`: Investigating the association between LPA LoF variants and APPOs. Run on an RStudio Workbench session.

## Running the Analysis (except UKB-RAP analysis)
From the main workspace, execute the analysis steps in order:
```bash
# Step 1: Get summary statistics
sbatch --account={XXXXX} scripts/step_1_get_sumstats.sh

# Step 2: Run MR analysis
sbatch --account={XXXXX} scripts/step_2_run_mr.sh

# Step 3: Sensitivity analyses
sbatch --account={XXXXX} scripts/step_3_other_sens.sh

# Step 4: IV checks
sbatch --account={XXXXX} scripts/step_4_iv_checks.sh
```
