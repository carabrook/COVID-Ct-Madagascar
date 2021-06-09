#!/bin/bash
#SBATCH --job-name=symp-gam
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

module load r/3.6.3
module load r-packages

Rscript pop-level-sym.R
