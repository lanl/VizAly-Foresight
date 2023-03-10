#!/bin/bash
### Visualization script ###

source $PROJWORK/evn_scripts/VizAly-CBench.bash.darwin
export PVPYTHON_PATH=/projects/exasky/ParaView-5.10.0-osmesa-MPI-Linux-Python3.9-x86_64/bin
export ORIGDATA_PATH=/projects/exasky/data/NYX/highz/512 #NVB_C009_l10n512_S12345T692_z42.hdf5
export DECOMPDATA_PATH=$PROJWORK/run_jinzhenw/decompressed

#img_orig=NVB_C009_l10n512_S12345T692_z54
#$PVPYTHON_PATH/pvpython vis_scripts/img_baryon_density_ReOr.py               $ORIGDATA_PATH/$img_orig.h5 $img_orig.h5 ${img_orig}_baryon_density.png 
#$PVPYTHON_PATH/pvpython vis_scripts/img_darkmatter_density_erdc_blue2gold.py $ORIGDATA_PATH/$img_orig.h5 $img_orig.h5 ${img_orig}_dark_matter_density.png 
#$PVPYTHON_PATH/pvpython vis_scripts/img_temperature_BuRd.py                  $ORIGDATA_PATH/$img_orig.h5 $img_orig.h5 ${img_orig}_temperature.png
#$PVPYTHON_PATH/pvpython vis_scripts/img_vx_cool2warm.py                      $ORIGDATA_PATH/$img_orig.h5 $img_orig.h5 ${img_orig}_velocity_x.png

export comp=mgard
for dname in $DECOMPDATA_PATH/$comp/*.h5;
do
dname=${dname##*/};
dname=${dname%.*};
#echo $dname;
#$PVPYTHON_PATH/pvpython vis_scripts/img_baryon_density_ReOr.py               $DECOMPDATA_PATH/$comp/$dname.h5 $dname.h5 syncDarwin/${comp}/relative_error/${dname}_baryon_density.png
#$PVPYTHON_PATH/pvpython vis_scripts/img_darkmatter_density_erdc_blue2gold.py $DECOMPDATA_PATH/$comp/$dname.h5 $dname.h5 syncDarwin/${comp}/relative_error/${dname}_dark_matter_density.png 
#$PVPYTHON_PATH/pvpython vis_scripts/img_temperature_BuRd.py    $DECOMPDATA_PATH/$comp/$dname.h5 $dname.h5  syncDarwin/$comp/relative_error/${dname}_temperature.png
$PVPYTHON_PATH/pvpython vis_scripts/img_vx_cool2warm.py        $DECOMPDATA_PATH/$comp/$dname.h5 $dname.h5  syncDarwin/$comp/relative_error/${dname}_velocity_x.png

done

#export comp=sz
#for dname in $DECOMPDATA_PATH/$comp/baryon_density/*.h5;
#do
#	dname=${dname##*/};
#	dname=${dname%.*};
#	echo $dname;
#	$PVPYTHON_PATH/pvpython vis_scripts/img_baryon_density_ReOr.py $DECOMPDATA_PATH/$comp/baryon_density/$dname.h5 $dname.h5 syncDarwin/${comp}/absolute_error/${dname}_baryon_density.png
#
#done
#
#for dname in $DECOMPDATA_PATH/$comp/dark_matter_density/*.h5;
#do
#	dname=${dname##*/};
#	dname=${dname%.*};
#	echo $dname;
#	$PVPYTHON_PATH/pvpython vis_scripts/img_darkmatter_density_erdc_blue2gold.py $DECOMPDATA_PATH/$comp/dark_matter_density/$dname.h5 $dname.h5 syncDarwin/${comp}/absolute_error/${dname}_dark_matter_density.png 
#
#done
#
#for dname in $DECOMPDATA_PATH/$comp/temperature/*.h5;
#do
#	dname=${dname##*/};
#	dname=${dname%.*};
#	echo $dname;
#        $PVPYTHON_PATH/pvpython vis_scripts/img_temperature_BuRd.py    $DECOMPDATA_PATH/$comp/temperature/$dname.h5 $dname.h5  syncDarwin/$comp/absolute_error/${dname}_temperature.png
#
#done
#
#for dname in $DECOMPDATA_PATH/$comp/velocity_x/*.h5;
#do
#	dname=${dname##*/};
#	dname=${dname%.*};
#	echo $dname;
#        $PVPYTHON_PATH/pvpython vis_scripts/img_vx_cool2warm.py        $DECOMPDATA_PATH/$comp/velocity_x/$dname.h5 $dname.h5  syncDarwin/$comp/absolute_error/${dname}_velocity_x.png
#
#done

