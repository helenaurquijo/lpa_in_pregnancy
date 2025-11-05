#!/bin/bash
#SBATCH --job-name=iv_checks
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --mem=200mb
#SBATCH --profile=task

# Exit immediately if any command fails
set -e

echo "Running script: $(basename $0)"
echo "Script run at: $(date '+%Y-%m-%d %H:%M:%S')"

source scripts/config.env
cd scripts

module load languages/R/4.3.3

## Run 14_pos_ctrl_MR.R
#  Performs colocalization analysis
echo "Running 14_pos_ctrl_MR.R..."
R CMD BATCH ${projectdir}/scripts/step_4_iv_checks/14_pos_ctrl_MR.R \
    ${projectdir}/logs/14_pos_ctrl_MR.Rout

# ... other IV checks ...