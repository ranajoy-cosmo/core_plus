#!/bin/bash -l
#SBATCH -p debug 
#SBATCH --qos premium
#SBATCH -N 20
#SBATCH -t 00:30:00
#SBATCH -J ts_sim
#SBATCH -o ts_sim_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
#export OMP_NUM_THREADS=1
export PYOPERATORS_NO_MPI=0

srun -n 120 ./sim_timestream.py config_sym4 run_mpi
