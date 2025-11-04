#!/bin/bash

# Script to recreate the lpa_in_pregnancy workspace structure
# Run this script from your desired project root directory

echo "Setting up lpa_in_pregnancy workspace..."

# Create main directory structure
mkdir -p archive
mkdir -p data/{UKB_RAP,exposures/{full_sumstats,instruments,raw},for_mrlap,outcomes/{full_sumstats,raw}}
mkdir -p logs
mkdir -p results/{UKB_RAP/{genetic_score,lof/plots},iv_checks,main_analysis,sensitivity_analyses/assoc_plots}
mkdir -p scripts/{UKB_RAP,functions,step_1_get_sumstats,step_2_run_mr,step_3_other_sens,step_4_iv_checks}

echo "Project structure:"
find . -type d | sort

