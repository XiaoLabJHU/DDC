#!/bin/bash
#SBATCH --job-name=matlab
#SBATCH --time=10:00:00
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
### Using more tasks because default memory is ~5GB per core
### 'shared' will share the node with other users
### 'parallel' use entire node (24,28,48, depends on node type)
### Try the script with --ntasks-per-node=1 and see what happens
#
#---------------------------------------------------------------------
# SLURM job script to run serial MATLAB
#---------------------------------------------------------------------

ml matlab
ml # confirm modules used
matlab -nodisplay -nosplash -nodesktop -r "Finish_Runs;"
echo "matlab exit code: $?"#!/bin/bash -l
