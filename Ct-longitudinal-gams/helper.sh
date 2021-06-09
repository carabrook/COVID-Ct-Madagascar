#!/bin/bash
#SBATCH --job-name=Ct-gams
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
module load r/3.6.3
module load r-packages
ht_helper.sh -m "r" -t RunScripts.sh -p 3
