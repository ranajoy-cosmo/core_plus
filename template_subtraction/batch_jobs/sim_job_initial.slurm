#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 45
#SBATCH -t 00:10:00
#SBATCH -J sim_init
#SBATCH -o out_files/sim_init.out

CONFIG_FILES_DIR=simulation.template_subtraction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.sim_config_a
srun -n 120 -N 15 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.sim_config_b
srun -n 120 -N 15 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.sim_config_TEMPLATE
srun -n 120 -N 15 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

wait
