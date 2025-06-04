#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o out_nbodysimpleNOGL/output.txt
#SBATCH -e out_nbodysimpleNOGL/errors.txt
./n_body $1 $2 $3