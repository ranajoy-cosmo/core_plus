#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 48 
#SBATCH -t 00:30:00
#SBATCH -J core_map_making
#SBATCH -o core_map_making_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
export PYOPERATORS_NO_MPI=0

srun -n 96 python map_maker_pol.py run_mpi
