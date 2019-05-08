#!/bin/bash -l

#SBATCH --qos=debug
#SBATCH -N 2
#SBATCH -c 1
#SBATCH --time=2
#SBATCH -A m2848
#SBATCH --constraint=haswell

projectPath=/project/projectdirs/m2848/lanl/VizAly-CBench

export HDF5_USE_FILE_LOCKING=FALSE

source $projectPath/scripts/VizAly-CBench.bash.cori
cd $projectPath/build/
srun -n 2 -c 1 $projectPath/build/CBench $projectPath/inputs/nyx_all.json