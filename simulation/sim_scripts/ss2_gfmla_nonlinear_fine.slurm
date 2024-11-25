#!/bin/bash

#SBATCH --mail-user=brichson@unc.edu  # email address
#SBATCH --mail-type=ALL  # Alerts sent when job begins, ends, or aborts
#SBATCH --nodes=1   # Number of nodes required
#SBATCH --ntasks=1   # Number of nodes required
#SBATCH --mem=16G  # Total memory (RAM) required, per node
#SBATCH --time=9-00:00:00  # Wall Clock time (dd-hh:mm:ss) [max of 14 days]
#SBATCH --array=0-9
#SBATCH --job-name=sim2_gfmla_nonlinear_fine # job name

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
module load r
Rscript ss2_gfmla_nonlinear_fine.R $SLURM_ARRAY_TASK_ID