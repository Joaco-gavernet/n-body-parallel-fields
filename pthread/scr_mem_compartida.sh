#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o out_memoria_compartida/output.txt
#SBATCH -e out_memoria_compartida/errors.txt
./n_body_pthreads $1 $2 $3 $4