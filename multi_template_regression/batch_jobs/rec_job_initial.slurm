#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 16
#SBATCH -t 00:04:00
#SBATCH -J rec_init
#SBATCH -o out_files/rec_init.out

CONFIG_FILES_DIR=simulation.multi_template_regression.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.rec_bolo_1a
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_bolo_1b
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_TEMPLATE_TD_1
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_TEMPLATE_SYNC_1
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_TEMPLATE_TD_1_QU
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_TEMPLATE_SYNC_1_QU
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_bolo_diff_QU_1
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_bolo_pair
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

wait
