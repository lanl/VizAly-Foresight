#!/bin/bash

for data in decompressed/*.hdf5;
do
	if [ ! -f  ${data%.*}.h5 ];
	then
		cp $data ${data%.*}.h5;
	fi
done

