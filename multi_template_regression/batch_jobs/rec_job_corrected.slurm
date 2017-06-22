#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 4
#SBATCH -t 00:10:00
#SBATCH -J rec_corr
#SBATCH -o out_files/rec_corr.out

CONFIG_FILES_DIR=simulation.multi_template_regression.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.rec_corrected
srun -n 96 -N 4 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

wait
