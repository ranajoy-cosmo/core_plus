#!/bin/bash -l
#SBATCH -p debug
#SBATCH --qos=premium
#SBATCH -N 40
#SBATCH -t 00:30:00
#SBATCH -J recon_map
#SBATCH -o recon_map_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
export PYOPERATORS_NO_MPI=0

srun -n 40 python map_maker_block.py config_mark_beam2 run_mpi
