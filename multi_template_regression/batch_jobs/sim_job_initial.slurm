#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 30
#SBATCH -t 00:10:00
#SBATCH -J sim_init
#SBATCH -o out_files/sim_init.out

CONFIG_FILES_DIR=simulation.multi_template_regression.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.sim_bolo_init
#srun -N 1 -n 1 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &
srun -N 15 -n 120 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.sim_TEMPLATE_init
#srun -N 1 -n 1 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &
srun -N 15 -n 120 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

wait
