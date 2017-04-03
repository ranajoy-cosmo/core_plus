#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 18
#SBATCH -t 00:12:00
#SBATCH -J sim_rescan
#SBATCH -o out_files/sim_rescan.out

CONFIG_FILES_DIR=simulation.pointing_correction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.resim_diff
srun -N 6 -n 48 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.resim_grad_co
srun -N 6 -n 48 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.resim_grad_cross
srun -N 6 -n 48 python $SIM_CODE_DIR/sim_timestream.py $CONFIG_FILE run_mpi &

wait
