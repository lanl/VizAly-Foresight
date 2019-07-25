The VizAly-CBench framework is based on the [factory design pattern](https://www.tutorialspoint.com/design_pattern/factory_pattern.htm). 
We expect that developers who contribute to this project will add new compressors and error metrics. 

### The metrics interface has the following functions:
* ``init`` to initialize the metrics interface
* ``execute`` executes the metrics
* ``close`` does cleanup

### Metrics Definition
The following error metrics being used are defined as:
  - absolute error
  	```
	  | actual_value - approx_value |
	```
  - relative error
  	```
	(| actual_value - approx_value |) / actual_value
	if  actual_value == 0
	   ...
	```
  - mean square error
  	```
	absolute error ^ 2
	```

  - minmax
  	```
	min_val = min(approx_value)
	max_val = max(approx_value)
	```
	
### Histogram
The following error metrics support a histogram output as a python script: absolute_error, relative_error, minmax. 
This allows for visualization of data distributions. Usage examples can be found in inputs/nyx/nyx_cbench_test.json

### Default Metrics
These metrics are always computed for all compressors:
  - Compression Ratio
    - compressed size / uncompressed size
  - Throughput
    - Size of data processed in MB / time take

	
