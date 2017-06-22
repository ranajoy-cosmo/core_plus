#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 64
#SBATCH -t 00:15:00
#SBATCH -J rec_aux
#SBATCH -o out_files/rec_aux.out

CONFIG_FILES_DIR=simulation.beam_correction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.rec_1a
srun -N 16 -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_1b
srun -N 16 -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_pair
srun -N 32 -n 96 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

wait
