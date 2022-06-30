#!/bin/bash
### Visualization script ###

source $PROJWORK/evn_scripts/VizAly-CBench.bash.darwin
export PVPYTHON_PATH=/projects/exasky/ParaView-5.10.0-osmesa-MPI-Linux-Python3.9-x86_64/bin
export ORIGDATA_PATH=/projects/exasky/data/NYX/highz/512 #NVB_C009_l10n512_S12345T692_z42.hdf5
export DECOMPDATA_PATH=$PROJWORK/run_jinzhenw/decompressed

for dname in $DECOMPDATA_PATH/sz_abs__*e-3__NVB_C009_l10n512_S12345T692_z54.h5;
do
	dname=${dname##*/};
	dname=${dname%.*};
	
	echo $dname;
	
	#$PVPYTHON_PATH/pvpython nyx_img_test_cold_hot.py $DECOMPDATA_PATH/$dname.h5 $dname.h5 ./figures/cold_hot/$dname.png
	$PVPYTHON_PATH/pvpython nyx_img_test.py $DECOMPDATA_PATH/$dname.h5 $dname.h5 ./figures/$dname.png


done


