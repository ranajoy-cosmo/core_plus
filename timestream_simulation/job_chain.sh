#!/bin/bash

sim=$(sbatch batch.slurm)
sim=${sim##* }

cd /global/homes/b/banerji/simulation/map_maker
sbatch -d afterok:$sim batch.slurm
