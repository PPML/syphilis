#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --array=1-5
#SBATCH --time=1-12:00:00
#SBATCH --mem-per-cpu=5GB
#SBATCH --job-name=optim
#SBATCH --error=sbatch-out/%x_%A_%a.err
#SBATCH --output=sbatch-out/%x_%A_%a.out

module load gcc/7.1.0-fasrc01
module load R/3.5.1-fasrc02
module load bzip2/1.0.6-fasrc01
srun Rscript optimize_with_overweighted_likelihoods.R ${SLURM_ARRAY_TASK_ID} $1 $2
