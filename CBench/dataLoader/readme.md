The VizAly-CBench framework is based on the [factory design pattern](https://www.tutorialspoint.com/design_pattern/factory_pattern.htm). 
We expect that developers who contribute to this project will add new compressors and error metrics. 

### The data loader interface has the following functions:
* ``init`` to initialize the compressor
* ``loadData`` to read in the data in parallel
* ``saveCompData`` to save out the compressed data as a binary file
* ``writeData`` to write out the decompressed data
* ``close`` to clean up
* ``satParam`` to pass parameters to the reafer
* ``loadUncompressedFields`` read in the scalars that do not get compressed



## Adding a new data loader
### Steps:
1. You first need to download, build and install the reader you want to add. Look at tools/externalDependencies/build_SZ.sh for an example script that does that. Add your script for your new compressor to that folder. 
2. Loaders should implement the dataLoaderInterface class located at CBench/dataLoader/dataLoaderInterface.hpp and have a CMakeLists.txt file. Look at the BLOSC example located at CBench/compressors/SZ.


```
4. Rerun cmake as follows ccamke ../CBench, turn USE_NEW_LOADER ON, a run configure. This should ask you to specify the locaton of your new reader's library and include path, e.g.:

5. Add your loader to the input json file and run the benchmark
