#!/bin/bash

#for data in decompressed/mgard/mgard*.hdf5;
#do
#	#if [ ! -f  ${data%.*}.h5 ];
#	#then
#	       cp $data ${data%.*}.h5;
#	#fi
#done

for data in decompressed/sz/velocity_x/*.hdf5;
do
	#if [ ! -f  ${data%.*}.h5 ];
	#then
	       cp $data ${data%.*}.h5;
	#fi
done

for data in decompressed/zfp/velocity_x/zfp*.hdf5;
do
	#if [ ! -f  ${data%.*}.h5 ];
	#then
	       cp $data ${data%.*}.h5;
	#fi
done
