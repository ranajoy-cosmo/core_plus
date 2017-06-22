#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 23
#SBATCH -t 00:06:00
#SBATCH -J sim_rescan
#SBATCH -o out_files/sim_rescan.out

CONFIG_FILES_DIR=simulation.multi_template_regression.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.sim_bolo_diff_1_rescan
srun -n 64 -N 8 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.sim_TEMPLATE_rescan
srun -n 120 -N 15 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

wait
