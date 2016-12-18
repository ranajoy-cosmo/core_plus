#!/bin/bash

sim_initial=$(sbatch sim_job_initial_2.slurm)
sim_initial=${sim_initial##* }

sbatch -d afterok:$sim_initial rec_job_initial_2.slurm
