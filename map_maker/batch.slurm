#!/bin/bash -l
#SBATCH -p debug 
#SBATCH -N 10
#SBATCH -t 00:10:00
#SBATCH -J recon_map
#SBATCH -o recon_map_%j.out

srun -n 240 python make_noise_mc.py 4800 
