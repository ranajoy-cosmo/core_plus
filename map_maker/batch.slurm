#!/bin/bash -l
#SBATCH -p regular
#SBATCH --qos=premium
#SBATCH -N 1
#SBATCH -t 00:02:00
#SBATCH -J recon_map
#SBATCH -o recon_map_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
export PYOPERATORS_NO_MPI=0

srun -n 1 python map_maker.py config_mark_beam1 run_mpi
