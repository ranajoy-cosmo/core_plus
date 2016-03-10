#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 26
#SBATCH -t 00:10:00
#SBATCH -J core_map_making
#SBATCH -o core_map_making_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
export PYOPERATORS_NO_MPI=0
#module switch PrgEnv-gnu/5.2.56 PrgEnv-intel

srun -n 104 -c 8 python map_maker_pol.py
