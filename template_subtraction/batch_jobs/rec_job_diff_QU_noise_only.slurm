#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 2
#SBATCH -t 00:10:00
#SBATCH -J rec_diff_QU_noise_only
#SBATCH -o out_files/rec_diff_QU_noise_only.out

CONFIG_FILES_DIR=simulation.template_subtraction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_diff_QU_noise_only

srun -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi
