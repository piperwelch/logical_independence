#!/bin/bash
# Specify a partition
#SBATCH --partition=bluemoon
# Request nodes 
#SBATCH --nodes=1
# Request some processor cores
#SBATCH --ntasks=1
# Maximum runtime
#SBATCH --time=05:00:00
# SBATCH --mem=16G                   # Total memory limit (e.g., 16 GB)
#SBATCH --ntasks=10              # Number of tasks (e.g., 4 tasks)
#SBATCH --cpus-per-task=2        # Number of CPU cores per task (2 cores per task)

python3 begin_evolving.py $1
