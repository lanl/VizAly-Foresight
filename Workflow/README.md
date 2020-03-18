# DRAW: Data Reduction Analysis Workflow tool

DRAW is a python tool for creating  Foresight workflows. It contains utilities for data reductio, analysis, and visualizing data



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
python3 -m <analysis_name> --input-file <absolute/relative path of input file>
```


For example,
```
# From Foresight root:
$ python3 -m Workflow.nyx.workflow --input-file ../inputs/nyx/darwin_nyx_galton_wflow_test.json
```


	

https://gaopinghuang0.github.io/2018/08/03/python3-import-and-project-layout
https://www.internalpointers.com/post/modules-and-packages-create-python-project