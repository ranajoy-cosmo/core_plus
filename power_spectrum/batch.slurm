#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 1 
#SBATCH -t 00:30:00
#SBATCH -J core_power_spectrum_estimation
#SBATCH -o core_power_spectrum_estimation_%j.out

export KMP_DUPLICATE_LIB_OK=TRUE
#module switch PrgEnv-gnu/5.2.56 PrgEnv-intel

srun -n 4 python make_power_spectra.py
