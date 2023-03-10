#!/bin/bash
## Running compression experiment

#salloc -N 1 -p skylake-gold

source $PROJWORK/evn_scripts/VizAly-CBench.bash.darwin
export PVPYTHON_PATH=/projects/exasky/ParaView-5.10.0-osmesa-MPI-Linux-Python3.9-x86_64/bin
export ORIGDATA_PATH=/projects/exasky/data/NYX/highz/512 #NVB_C009_l10n512_S12345T692_z42.hdf5
export DECOMPDATA_PATH=$PROJWORK/run_jinzhenw/decompressed

mpirun -np 4 ../build/CBench ../inputs/nyx/nyx_img_compression_sz_abs_baryon_density_test.json
mpirun -np 4 ../build/CBench ../inputs/nyx/nyx_img_compression_zfp_abs_baryon_density_test.json
mpirun -np 4 ../build/CBench ../inputs/nyx/nyx_img_compression_sz_abs_temperature_test.json
mpirun -np 4 ../build/CBench ../inputs/nyx/nyx_img_compression_zfp_abs_temperature_test.json

