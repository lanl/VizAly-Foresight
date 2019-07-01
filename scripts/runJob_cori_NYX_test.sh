#!/bin/bash -l

#SBATCH --qos=debug
#SBATCH -N 2
#SBATCH -c 1
#SBATCH --time=2
#SBATCH -A m2294
#SBATCH --constraint=haswell

projectPath=$(pwd)

export HDF5_USE_FILE_LOCKING=FALSE
source $projectPath/scripts/VizAly-CBench.bash.cori

cd $projectPath/build/
srun -n 2 -c 1 $projectPath/build/CBench $projectPath/inputs/nyx/nyx_cbench_test.json