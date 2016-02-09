#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 48
#SBATCH -t 00:30:00
#SBATCH -J core_ts_sim
#SBATCH -o core_ts_sim_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
#module switch PrgEnv-gnu/5.2.56 PrgEnv-intel
export OMP_NUM_THREADS=1

srun -n 192 -c 6 ./sim_timestream_pol.py
