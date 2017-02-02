#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 8
#SBATCH -t 00:10:00
#SBATCH -J rec_corr
#SBATCH -o out_files/rec_corr.out

CONFIG_FILES_DIR=simulation.template_subtraction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_corrected_1
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_corrected_2
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_corrected
srun -n 96 -N 4 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

wait
