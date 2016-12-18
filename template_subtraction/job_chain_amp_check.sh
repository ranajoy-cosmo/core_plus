#!/bin/bash

BATCH_DIR=$HOME/simulation/template_subtraction/batch_jobs

sim_initial=$(sbatch $BATCH_DIR/sim_job_initial.slurm)
sim_initial=${sim_initial##* }

rec_initial=$(sbatch -d afterok:$sim_initial $BATCH_DIR/rec_job_amp_check.slurm)
rec_initial=${rec_initial##* }
