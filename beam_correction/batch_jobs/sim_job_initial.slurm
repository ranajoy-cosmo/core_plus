#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 112
#SBATCH -t 00:30:00
#SBATCH -J sim_init
#SBATCH -o out_files/sim_init.out

CONFIG_FILES_DIR=simulation.beam_correction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.sim_1
srun -N 64 -n 384 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.sim_tm_grad
srun -N 48 -n 384 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

wait
