#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=5GB
#SBATCH --job-name=run_interventions
#SBATCH --error=sbatch-out/%x_%A_%a.err
#SBATCH --output=sbatch-out/%x_%A_%a.out

module load gcc/7.1.0-fasrc01
module load R/3.5.1-fasrc02
module load bzip2/1.0.6-fasrc01
srun Rscript run_interventions.R $1
