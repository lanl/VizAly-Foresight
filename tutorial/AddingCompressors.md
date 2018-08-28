In this tutorial, we will go over how to a compressor to CBench. The example we will use in this case is zfp - https://github.com/LLNL/zfp

## Steps:
1. You first need to download, build and install zfp. Use the build_zfp.sh script for that. 
2. Compressors should implement the CompressorInterface class located at src/compressors/compressorInterface.hpp
    1. Add the zfp folder containing zfpCompressor.hpp, a class implementing the CompressorInterface for zfp and a CMakeLists.txt file to the compressors folder
3. Modify main.cpp in src by adding contents from main_modifications.diff
4. Rerun cmake as follows ccamke ../src, turn USE_ZFP ON, a run configure. This should ask you to specify the locaton of zfp library and include path, e.g.:
    1. ZFP_LIBRARY     : ../ExternalDependencies/zfp/install/lib/libzfp.so 
    2. ZFP_INCLUDE_PATH: ../ExternalDependencies/zfp/install/include
    3. Configure, Generate, and rerun make
5. Run with the zfp.json file as input e.g. mpirun -np 2 CBench ../tutorial/zfp_tutorial_files/zfp.json

