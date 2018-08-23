In this tutorial, we will go over how to a compressor to CBench. The example we will use in this case is zfp - https://github.com/LLNL/zfp

## Steps:
1. First download, build and install zfp. Once installed, the install path for zfp should have a lib and an include folder.
2. Compressors should implement the CompressorInterface class located at src/compressors/compressorInterface.hpp
    1. Add zfpCompressor.hpp, a class implementing the CompressorInterface for zfp to the compressor folder
3. Modify the CMakeLists.txt file in src to allow it to link to zfp
    1. Add the contents of zfp_CMakeListsAddition.txt to the CMakeLists.txt file in src beofre the find_package(MPI) line 
4. Modify main.cpp in src by adding contents from main_modifications.diff
5. Rerun cmake as follows ccamke ../src, turn USE_ZFP ON, a run configure. This should ask you to specify the locaton of zfp library and include path, e.g.:
    1. ZFP_LIBRARY     : /home/pascal/software/zfp-0.5.3/install/lib/libzfp.so 
    2. ZFP_INCLUDE_PATH: /home/pascal/software/zfp-0.5.3/install/include
    3. Configure, Generate, and rerun make
6. Run with the zfp.json file as input e.g. mpirun -np 2 CBench ../tutorial/zfp_tutorial_files/zfp.json

