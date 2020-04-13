#!/bin/bash

# VTK
echo "Building VTK ... "

echo $PLATFORM

if [ "$PLATFORM" = "travis" ]; then
	echo "Travis: Using docker VTK build"
else
	git clone https://gitlab.kitware.com/vtk/vtk.git -b v8.2.0
	cd vtk
	mkdir install
	mkdir build
	cd build
	cmake .. -DVTK_Group_MPI:BOOL=true -DBUILD_TESTING:BOOL=false -DCMAKE_INSTALL_PREFIX=../install
	make -j
	make install
	cd ..
	cd ..
fi

echo "Building VTK done!"
