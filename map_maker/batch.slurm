#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 128
#SBATCH -t 00:10:00
#SBATCH -J core_map_making
#SBATCH -o core_map_making_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
export PYOPERATORS_NO_MPI=0

srun -n 1024 python map_maker_pol.py
