#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 60
#SBATCH -t 00:30:00
#SBATCH -J core_ts_sim
#SBATCH -o core_ts_sim_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
#module switch PrgEnv-gnu/5.2.56 PrgEnv-intel
export OMP_NUM_THREADS=8
export PYOPERATORS_NO_MPI=0

srun -n 180 -c 8 ./sim_timestream_pol.py
