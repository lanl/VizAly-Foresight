# PAT: Python Analysis Toolkit

This folder contains PAT, the Python package that is used to run the Foresight workflows. The Analysis toolkit will generate a Cinema database from the input file specified.


## Prerequisites:

Minimum Requirements
* Python 3.6 or higher with (+matplotlib=3.0.2, +apsw=3.9.2, +numpy=1.15.4)
* SLURM (for job launching)


## Running Foresight

Foresight can be run in 3 different modes:
1. Full workflow (Run compression benchmark, analysis on decompressed data, and Visualize the results). This is the default behavior.
2. Analysis (Analysis on decompressed data, and Visualize the results)
3. Cinema (Visualize the results)


### Running Full analysis
This will run the full analysis workflow and generate a cinema database. 

To generate the analysis for nyx, use the command as follows:
```
python3 -m <analysis_name> --input-file <absolute path of input file>
```


For example,
```
cd Analysis
$ python3 -m pat.nyx.workflow --input-file ../inputs/nyx/darwin_nyx_galton_wflow_test.json
```


### Options
Append the following options for the behavior below
* --preview: only create the scripts but do not launch them

* --analysis_cinema: opmits CBench, existing output files must be provided
	* metrics file must be proivided in ["inputs"] e.g.
		```
		"input": 
		{
			"metrics-csv" : "/projects/VizAly-Foresight/test/metrics_.csv",
			.
			.
			.
		}
		```
	* paths of files must be provided ["pat"] e.g.
		```
		"pat" :
		{
			.
			.
			.

			"input-files": [
	            {
	                "output-prefix": "orig",
	                "path": "/projects/VizAly-Foresight/testing/data/z255_32.h5"
	            },
	            {
	                "output-prefix": "SZ_",
	                "path": "/projects/VizAly-Foresight/test/cbench/decompressed_files/SZ___z255_32.h5"
	            }
	        ]
	    }
	    ```
* --cinema: only produces cinema output from existing CBench and analysis runs. 
	* metrics and analysis-results location must be proivided in ["inputs"] and  e.g.
		```
			"input": 
			{
				"metrics-csv" : "/projects/VizAly-Foresight/test/metrics_.csv",
				"analysis-results" : "/projects/VizAly-Foresight/testB/",
				.
				.
				.
			}
		```
	* paths of files must be provided e.g.
		```
		"pat" :
		{
			.
			.
			.

			"input-files": [
	            {
	                "output-prefix": "orig",
	                "path": "/projects/VizAly-Foresight/testing/data/z255_32.h5"
	            },
	            {
	                "output-prefix": "SZ_",
	                "path": "/projects/VizAly-Foresight/test/cbench/decompressed_files/SZ___z255_32.h5"
	            }
	        ]
	    }
	    ```