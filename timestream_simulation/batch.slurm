#!/bin/bash -l
#SBATCH -p regular 
#SBATCH --qos premium
#SBATCH -N 3
#SBATCH -t 00:30:00
#SBATCH -J ts_sim
#SBATCH -o ts_sim_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
#export OMP_NUM_THREADS=1
export PYOPERATORS_NO_MPI=0

srun -n 24 ./sim_timestream.py config_350_T run_mpi
