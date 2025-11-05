#!/bin/bash

# Script to recreate the lpa_in_pregnancy workspace structure
# Run this script from your desired project root directory

echo "Setting up lpa_in_pregnancy workspace..."

# Create main directory structure
mkdir -p ../data/{UKB_RAP,exposures/{full_sumstats,instruments,raw},for_mrlap,outcomes/{full_sumstats,raw}}
mkdir -p ../logs
mkdir -p ../results/{UKB_RAP/{genetic_score,lof/plots},iv_checks,main_analysis,sensitivity_analyses/assoc_plots}

echo "Project structure:"
find . -type d | sort