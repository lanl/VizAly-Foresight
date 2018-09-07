The VizAly-CBench framework is based on the [factory design pattern](https://www.tutorialspoint.com/design_pattern/factory_pattern.htm). 
We expect that developers who contribute to this project will add new compressors and error metrics. 

### The metrics interface has the following functions:
* ``init`` to initialize the metrics interface
* ``execute`` executes the metrics
* ``close`` does cleanup


The following error metrics being used are defined as:
  - absolute error
	- | actual_value - approx_value |
  - relative error
	- (| actual_value - approx_value |) / actual_value
	- if  actual_value == 0
	  - ...
  - mean square error
	- absolute error ^ 2


Static Metric computed for all compressors:
  - Compression Ration
    - compressed size / uncompressed size
  - Throughput
    - Size of data processed in MB / time take

	