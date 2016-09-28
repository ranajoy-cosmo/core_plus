#!/bin/bash -l
#SBATCH -p debug
#SBATCH --qos=premium
#SBATCH -N 3
#SBATCH -t 00:15:00
#SBATCH -J recon_map
#SBATCH -o recon_map_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
export PYOPERATORS_NO_MPI=0

srun -n 24 python map_maker.py config_recon_test run_mpi
