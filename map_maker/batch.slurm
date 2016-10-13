#!/bin/bash -l
#SBATCH -p regular
#SBATCH --qos=premium
#SBATCH -N 5
#SBATCH -t 00:15:00
#SBATCH -J recon_map
#SBATCH -o recon_map_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
export PYOPERATORS_NO_MPI=0

srun -n 40 python map_maker_signal_diff.py config_signal run_mpi
