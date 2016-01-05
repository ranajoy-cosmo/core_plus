#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 3
#SBATCH -t 00:30:00
#SBATCH -J core_map_making
#SBATCH -o core_map_making_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
module switch PrgEnv-gnu/5.2.56 PrgEnv-intel

srun -n 36 -c 2 python map_maker_pol.py
