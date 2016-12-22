#!/bin/bash -l
#SBATCH -p debug
#SBATCH --qos=premium
#SBATCH -N 11
#SBATCH -t 00:20:00
#SBATCH -J recon_map
#SBATCH -o recon_map_%j.out

srun -n 32 python map_maker.py config1 run_mpi
