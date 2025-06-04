#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH -o out_hybrid/output.txt
#SBATCH -e out_hybrid/errors.txt
mpirun --bind-to none ./hybrid $1 $2 $3 $4