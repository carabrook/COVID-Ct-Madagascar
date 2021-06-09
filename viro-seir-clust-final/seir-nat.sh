#!/bin/bash
#SBATCH --job-name=cluseirNat
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00

singularity run rlsoda_0.1.sif Rscript fit_seir_national.R
