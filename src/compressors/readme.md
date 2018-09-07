The VizAly-CBench framework is based on the [factory design pattern](https://www.tutorialspoint.com/design_pattern/factory_pattern.htm). 
We expect that developers who contribute to this project will add new compressors and error metrics. 

### The compressor interface has the following functions:
* ``init`` to initialize the compressor
* ``compress`` to compress the data
* ``decompress`` to decompress the data
* ``close`` to remove anything after the test



## Adding a new compressor
### Steps:
1. You first need to download, build and install the compressor you want to add. Look at tools/externalDependencies/build_SZ.sh for an example script that does that. Add your script for your new compressor to that folder. 
2. Compressors should implement the CompressorInterface class located at src/compressors/compressorInterface.hpp and have a CMakeLists.txt file. Look at the BLOSC example located at src/compressors/SZ.
3. Modify src/compressors/compressorFactory.hpp to add your new compressor. 
 ``
	#ifdef NEW_COMPRESSOR
	
        else if (compressorName == "NewCompressor")
	
          return new NewCompressor();
	  
      #endif
 ``
4. Rerun cmake as follows ccamke ../src, turn USE_NEW_COMPRESSOR ON, a run configure. This should ask you to specify the locaton of your new compressor's library and include path, e.g.:
    1. NEW_COMPRESSOR_LIBRARY     : ../<path here>/libName.so 
    2. NEW_COMPRESSOR_INCLUDE_PATH: ../<path here>/include
    3. Configure, Generate, and rerun make
5. Add your compressor to the input json file and run the benchmark
