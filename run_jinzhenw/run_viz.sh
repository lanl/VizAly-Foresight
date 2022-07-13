#!/bin/bash
### Visualization script ###

source $PROJWORK/evn_scripts/VizAly-CBench.bash.darwin
export PVPYTHON_PATH=/projects/exasky/ParaView-5.10.0-osmesa-MPI-Linux-Python3.9-x86_64/bin
export ORIGDATA_PATH=/projects/exasky/data/NYX/highz/512 #NVB_C009_l10n512_S12345T692_z42.hdf5
export DECOMPDATA_PATH=$PROJWORK/run_jinzhenw/decompressed

for dname in $DECOMPDATA_PATH/sz_rel__*-4*.h5;
do
	dname=${dname##*/};
	dname=${dname%.*};
	
	echo $dname;
	
	#$PVPYTHON_PATH/pvpython z54_byrondensity_cool2hot_singlecolor.py $DECOMPDATA_PATH/$dname.h5 $dname.h5 ${dname}_baryon_density.png 
        $PVPYTHON_PATH/pvpython z54_temp_cool2warm_threecolor.py         $DECOMPDATA_PATH/$dname.h5 $dname.h5 ${dname}_temperature.png
        $PVPYTHON_PATH/pvpython z54_vx_cool2warm_threecolor.py           $DECOMPDATA_PATH/$dname.h5 $dname.h5 ${dname}_velocity_x.png



done


