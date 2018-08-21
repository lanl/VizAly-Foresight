In this tutorial, we will go over how to a compressor to CBench. The example we will use in this case is zfp - https://github.com/LLNL/zfp

##Steps:
1. First download, build and install zfp. Once installed, the install path for zfp should have a lib and an include folder.
2. Compressors shoudl implement the CompressorInterface class located at compressors/compressorInterface.hpp
	2.1 Add a subfolder zfp under src/compressors
	2.1 An example of that for zfp is available in the zfp folder in the zfp_tutorial_files
3. Modify the CMakeLists.txt file in src to allow it to link to zfp
	3.1 Add the contents of zfp_CMakeListsAddition.txt to the CMakeLists.txt file in src beofre the find_package(MPI) line 
4. Modify main.cpp in src by adding contents from main_modifications.diff
5. Rerun cmake as follows ccamke ../src, turn USE_ZFP ON, a run configure. This shoudl ask you to specify the locaton of zfp library and include path
	ZFP_LIBRARY should point to libzfp.so e.g ZFP_LIBRARY /home/pascal/software/zfp-0.5.3/install/lib/libzfp.so 
	ZFP_INCLUDE_PATH should point to the include folders e.g. ZFP_INCLUDE_PATH: /home/pascal/software/zfp-0.5.3/install/include
  Configure and Generate and rerun make as follows: make -j
6. Run with the zfp.json file as input e.g. mpirun -np 2 CBench ../tutorial/zfp_tutorial_files/zfp.json

