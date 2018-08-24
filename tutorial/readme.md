The VizAly-CBench framework is based on the [factory design pattern](https://www.tutorialspoint.com/design_pattern/factory_pattern.htm). 
We expect that developers who contribute to this project will add new compressors and error metrics. 

### The compressor interface has the following functions:
* ``init`` to initialize the compressor
* ``compress`` to compress the data
* ``decompress`` to decompress the data
* ``close`` to remove anything after the test

AddingCompressors.md shows how to add a compressor to this framework.

### The metrics interface has the following functions:
* ``init`` to initialize the metrics interface
* ``execute`` executes the metrics
* ``close`` does cleanup
