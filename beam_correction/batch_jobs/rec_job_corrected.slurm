#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 12
#SBATCH -t 00:15:00
#SBATCH -J rec_corr
#SBATCH -o out_files/rec_corr.out

CONFIG_FILES_DIR=simulation.beam_correction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.rec_corrected
srun -N 12 -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi
