#!/bin/bash

#SBATCH --job-name=mkSNO
#SBATCH --partition=batch
#SBATCH --output=mksno-slurm-%N.%A.%a.out
#SBATCH --error=mksno-slurm-%N.%A.%a.err
#SBATCH --qos=normal
#SBATCH --account=pi_strow
#SBATCH -N1
#SBATCH --mem-per-cpu=10000
#SBATCH --cpus-per-task 2
#SBATCH --array=1

MATLAB='/usr/cluster/matlab/2016b/bin/matlab'
MATOPTS=' -nodisplay -nojvm -nosplash'

srun $MATLAB $MATOPTS -r "run_makeSNO_batch; exit"
