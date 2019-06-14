#!/bin/bash
#SBATCH -N 8
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH --mail-user=pascalgrosset@lanl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 00:20:00
#SBATCH -A m2848


# Set environment
projectPath=/project/projectdirs/m2848/lanl/VizAly-CBench

export HDF5_USE_FILE_LOCKING=FALSE
source $projectPath/scripts/VizAly-CBench.bash.cori


#run the application:
cd $projectPath/build/
srun -n 64 -c 8 $projectPath/build/CBench $projectPath/inputs/nyx/cori_nyx_thermal_2048_all.json