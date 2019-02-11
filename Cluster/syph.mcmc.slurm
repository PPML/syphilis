#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-11:00
#SBATCH -p serial_requeue
#SBATCH --mem 600
#SBATCH -J R_test
#SBATCH -o R_test.out
#SBATCH -e R_test.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=atuite@hsph.harvard.edu

R CMD BATCH --quiet --no-save --no-restore syph_mcmc.R mcmc.out
