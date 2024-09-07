#!/bin/bash

#SBATCH --account=rubin:developers      # Account name
#SBATCH --partition=milano
#
#SBATCH --job-name=four_cycle
#SBATCH --output=output-%j.txt
#SBATCH --error=output-%j.txt
#
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20g
#
#SBATCH --time=1-14:10:00
#SBATCH --chdir=/sdf/data/rubin/shared/fbs_sims/sims_featureScheduler_runs3.6/four_cycle
#SBATCH --mail-user=yoachim@uw.edu
#
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yoachim@uw.edu

## Set up the evironment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate rubin
export OPENBLAS_NUM_THREADS=1

python four_cycle.py 
