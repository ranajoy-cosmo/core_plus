#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH -J y_est
#SBATCH -o out_files/y_est.out

CONFIG_FILES_DIR=simulation.template_subtraction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.default_config

srun -n 24 python $REC_CODE_DIR/y_estimator.py $CONFIG_FILE 120
