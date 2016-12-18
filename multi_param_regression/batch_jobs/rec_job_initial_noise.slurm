#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 20
#SBATCH -t 00:10:00
#SBATCH -J rec_init
#SBATCH -o out_files/rec_init.out

CONFIG_FILES_DIR=simulation.template_subtraction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_a
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_a_noise
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_b
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_b_noise
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_TEMPLATE
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_diff_QU
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_diff_QU_noise
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_TEMPLATE_QU
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_pair
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_pair_noise
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

wait
