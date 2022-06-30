#!/bin/bash

for data in decompressed/sz_abs__*e-3__NVB_C009_l10n512_S12345T692_z54.hdf5;
do
	#if [ ! -f  ${data%.*}.h5 ];
	#then
	       cp $data ${data%.*}.h5;
	#fi
done

