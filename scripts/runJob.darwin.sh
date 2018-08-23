#!/bin/bash
#SBATCH -N 8
#SBATCH --ntasks-per-node 8
#SBATCH -p scaling

# load modules
source /projects/exasky/VizAly-CBench/scripts/VizAly-CBench.bash.darwin

# go to folder 
cd /projects/exasky/VizAly-CBench/build

# Run:
mpirun /projects/exasky/VizAly-CBench/build/CBench /projects/exasky/VizAly-CBench/inputs/all_step499.json

