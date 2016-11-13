#!/bin/bash

sim_a=$(qsub sim_job_a.pbs)
sim_b=$(qsub sim_job_b.pbs)
sim_TEMPLATE=$(qsub sim_job_TEMPLATE.pbs)

rec_diff=$(qsub -W depend=afterok:$sim_a:$sim_b rec_job_diff.pbs)
rec_TEMPLATE_QU=$(qsub -W depend=afterok:$sim_TEMPLATE rec_job_TEMPLATE_QU.pbs)

sim_rescan_diff=$(qsub -W depend=afterok:$rec_diff sim_job_rescan_diff.pbs)
sim_rescan_TEMPLATE=$(qsub -W depend=afterok:$rec_TEMPLATE_QU sim_job_rescan_TEMPLATE.pbs)

TEMPLATE_amp_estimate=$(qsub -W depend=afterok:$sim_rescan_TEMPLATE:$sim_rescan_diff job_TEMPLATE_estimation.pbs)

rec_corrected=$(qsub -W depend=afterok:$TEMPLATE_amp_estimate rec_job_corrected.pbs)

qsub -W depend=afterok:$sim_a rec_job_a.pbs
qsub -W depend=afterok:$sim_b rec_job_b.pbs
qsub -W depend=afterok:$sim_TEMPLATE rec_job_TEMPLATE.pbs
qsub -W depend=afterok:$sim_a:$sim_b rec_job_pair.pbs
