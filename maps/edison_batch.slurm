#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -J core_gen_maps
#SBATCH -o core_gen_maps_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
#module switch PrgEnv-gnu/5.2.56 PrgEnv-intel
export OMP_NUM_THREADS=1

srun -n 1 ./generate_map_from_spectra.py
