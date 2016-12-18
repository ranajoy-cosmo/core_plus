#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 14
#SBATCH -t 00:10:00
#SBATCH -J rec_init
#SBATCH -o out_files/rec_init.out

CONFIG_FILES_DIR=simulation.template_subtraction.config_files
SIM_CODE_DIR=$HOME/simulation/timestream_simulation
REC_CODE_DIR=$HOME/simulation/map_maker

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_pair_1
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_pair_2
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

#CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_pair_3
#srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_all
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_1a
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_1b
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_2a
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_2b
srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

#CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_3a
#srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &
#
#CONFIG_FILE=$CONFIG_FILES_DIR.rec_config_3b
#srun -n 48 -N 2 python $REC_CODE_DIR/map_maker.py $CONFIG_FILE run_mpi &

wait
