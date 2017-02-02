#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 60
#SBATCH -t 00:06:00
#SBATCH -J sim_rescan
#SBATCH -o out_files/sim_rescan.out

CONFIG_FILES_DIR=simulation.template_subtraction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.sim_config_rescan_diff_1
srun -n 120 -N 15 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.sim_config_rescan_TEMPLATE_1
srun -n 120 -N 15 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.sim_config_rescan_diff_2
srun -n 120 -N 15 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.sim_config_rescan_TEMPLATE_2
srun -n 120 -N 15 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

wait
