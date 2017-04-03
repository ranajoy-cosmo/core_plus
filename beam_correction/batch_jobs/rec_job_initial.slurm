#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 48
#SBATCH -t 00:10:00
#SBATCH -J rec_init
#SBATCH -o out_files/rec_init.out

CONFIG_FILES_DIR=simulation.beam_correction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.rec_diff
srun -N 12 -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &
#
#CONFIG_FILE=$CONFIG_FILES_DIR.rec_tm_grad_co
#srun -N 12 -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &
#
#CONFIG_FILE=$CONFIG_FILES_DIR.rec_tm_grad_cross
#srun -N 12 -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &
#
CONFIG_FILE=$CONFIG_FILES_DIR.rec_tm_grad_co_co
srun -N 12 -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_tm_grad_cross_cross
srun -N 12 -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_tm_grad_co_cross
srun -N 12 -n 48 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

wait
