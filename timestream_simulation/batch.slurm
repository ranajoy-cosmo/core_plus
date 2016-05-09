#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 1 
#SBATCH -t 00:30:00
#SBATCH -J core_ts_sim
#SBATCH -o core_ts_sim_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
#export OMP_NUM_THREADS=1
export PYOPERATORS_NO_MPI=0

srun -n 8 ./sim_timestream_pol.py run_mpi
