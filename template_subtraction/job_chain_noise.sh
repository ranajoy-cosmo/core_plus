#!/bin/bash

BATCH_DIR=$HOME/simulation/template_subtraction/batch_jobs

sim_initial=$(sbatch $BATCH_DIR/sim_job_initial.slurm)
sim_initial=${sim_initial##* }

rec_initial=$(sbatch -d afterok:$sim_initial $BATCH_DIR/rec_job_initial_noise.slurm)
rec_initial=${rec_initial##* }

sim_rescan=$(sbatch -d afterok:$rec_initial $BATCH_DIR/sim_job_rescan.slurm)
sim_rescan=${sim_rescan##* }

y_est=$(sbatch -d afterok:$sim_rescan $BATCH_DIR/y_estimation.slurm)
y_est=${y_est##* }

sbatch -d afterok:$y_est $BATCH_DIR/rec_job_corrected.slurm 
