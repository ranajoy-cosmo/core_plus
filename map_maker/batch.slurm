#!/bin/bash -l
#SBATCH -p debug 
#SBATCH -N 22
#SBATCH -t 00:15:00
#SBATCH -J recon_map
#SBATCH -o recon_map_%j.out

CONFIG_FILE_FOLDER=simulation.map_maker.config_files

CONFIG_FILE=$CONFIG_FILE_FOLDER.config1a
srun -N 22 -n 64 python map_maker.py $CONFIG_FILE run_mpi
#srun -N 22 -n 64 python map_maker.py $CONFIG_FILE run_mpi &
#
#CONFIG_FILE=$CONFIG_FILE_FOLDER.config1b
#srun -N 22 -n 64 python map_maker.py $CONFIG_FILE run_mpi &
#
#CONFIG_FILE=$CONFIG_FILE_FOLDER.config1a1b
#srun -N 44 -n 128 python map_maker.py $CONFIG_FILE run_mpi &
#
#wait
