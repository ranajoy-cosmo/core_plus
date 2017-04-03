#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 24
#SBATCH -t 00:30:00
#SBATCH -J recon_map
#SBATCH -o recon_map_%j.out

CONFIG_FILE_FOLDER=simulation.map_maker.config_files

CONFIG_FILE=$CONFIG_FILE_FOLDER.config_test
srun -n 192 python map_maker.py $CONFIG_FILE run_mpi
