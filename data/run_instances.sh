#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=60:00:00
#SBATCH --array=1-1
#SBATCH --output=../data/%A_%a.out
#SBATCH --account=def-vidalthi

module load gurobi/12.0.0
module load python

virtualenv --no-download  ~/env_gurobi
source ~/env_gurobi/bin/activate
python read_instance.py