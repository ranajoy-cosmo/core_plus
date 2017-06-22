#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:15:00
#SBATCH -J y_est
#SBATCH -o out_files/y_est.out

CONFIG_FILES_DIR=simulation.beam_correction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.y_est_config
srun -n 24 python multi_y_estimator.py $CONFIG_FILE 768
