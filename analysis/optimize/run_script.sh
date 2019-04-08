#!/bin/bash
#SBATCH -p shared
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --array=1-25
#SBATCH --time=1-12:00:00
#SBATCH --mem-per-cpu=5GB
#SBATCH --job-name=optim
#SBATCH --error=sbatch-out/%x_%A_%a.err
#SBATCH --output=sbatch-out/%x_%A_%a.out

source ~/load-gcc-and-R.sh
srun Rscript optimize.R $1 ${SLURM_ARRAY_TASK_ID}
