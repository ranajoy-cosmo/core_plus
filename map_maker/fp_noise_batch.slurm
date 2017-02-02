#!/bin/bash -l
#SBATCH -p debug 
#SBATCH -N 40
#SBATCH -t 00:30:00
#SBATCH -J recon_map
#SBATCH -o recon_map_%j.out

CONFIG_FILE_FOLDER=simulation.map_maker.config_files

CONFIG_FILE=$CONFIG_FILE_FOLDER.config_fp_noise
srun -n 480 python map_maker.py $CONFIG_FILE run_mpi
