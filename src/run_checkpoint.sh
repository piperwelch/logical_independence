#!/bin/bash
# Specify a partition
#SBATCH --partition=bluemoon
# Request nodes 
#SBATCH --nodes=1
# Request some processor cores
#SBATCH --ntasks=1
# Maximum runtime
#SBATCH --time=30:00:00

python3 checkpoint.py $1 $2