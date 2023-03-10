#!/bin/bash
#
#SBATCH --job-name=sz_zfp_velocity_x
#SBATCH --output=sz_zfp_velocity_x.log
#SBATCH --partition=skylake-gold
#SBATCH --time=60:00
#SBATCH --nodes=5

source $PROJWORK/evn_scripts/VizAly-CBench.bash.darwin
export PVPYTHON_PATH=/projects/exasky/ParaView-5.10.0-osmesa-MPI-Linux-Python3.9-x86_64/bin
export ORIGDATA_PATH=/projects/exasky/data/NYX/highz/512 #NVB_C009_l10n512_S12345T692_z42.hdf5
export DECOMPDATA_PATH=$PROJWORK/run_jinzhenw/decompressed

#mpirun -np 20 ../build/CBench ../inputs/nyx/nyx_img_compression_sz_abs_baryon_density.json
#mpirun -np 20 ../build/CBench ../inputs/nyx/nyx_img_compression_sz_abs_dark_matter_density.json
#mpirun -np 20 ../build/CBench ../inputs/nyx/nyx_img_compression_zfp_abs_baryon_density.json
#mpirun -np 20 ../build/CBench ../inputs/nyx/nyx_img_compression_zfp_abs_dark_matter_density.json
#mpirun -np 20 ../build/CBench ../inputs/nyx/nyx_img_compression_sz_abs_temperature.json
mpirun -np 20 ../build/CBench ../inputs/nyx/nyx_img_compression_sz_abs_velocity_x.json
#mpirun -np 20 ../build/CBench ../inputs/nyx/nyx_img_compression_zfp_abs_temperature.json
mpirun -np 20 ../build/CBench ../inputs/nyx/nyx_img_compression_zfp_abs_velocity_x.json
#mpirun -np 4 ../build/CBench ../inputs/nyx/nyx_img_compression_zfp_abs_jw.json
#mpirun -np 20 ../build/CBench ../inputs/nyx/nyx_img_compression_MGARD.json

#mpirun -np 1 ./run_viz.sh
