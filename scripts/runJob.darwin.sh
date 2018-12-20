#!/bin/bash
#SBATCH -N 8
#SBATCH --ntasks-per-node 8
#SBATCH -p scaling

projectPath=/projects/exasky/VizAly-CBench/

# load modules
source $projectPath/scripts/VizAly-CBench.bash.darwin

# go to folder 
cd $projectPath/build

# Run:
mpirun $projectPath/build/CBench $projectPath/inputs/darwin_step499.json

