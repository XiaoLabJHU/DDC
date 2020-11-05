#!/bin/bash -l

#SBATCH
#SBATCH --job-name=long_time
#SBATCH --time=50:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-type=end
#SBATCH --mail-user=cbohrer1@jhu.edu

module load matlab     #### load matlab module, by defult 2015

matlab -nodesktop -nosplash -r Finish_Runs > long2.log