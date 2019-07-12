# PAT: Python Analysis Toolkit

This folder contains PAT, the Python package that is used to run the Foresight workflows. The Analysis toolkit will generate a Cinema database from the input file specified.

## Prerequisites:

* Python 2.7

There is a `requirements.txt` file in the top-level of this repository. All of the required Python packages can be installed as: `pip install -r requirements.txt`

## Running Foresight
This will run the full analysis workflow and generate a cinema database. 

To run the anlysis for nyx, the command is as follows:
```
python <analysis_name> --input-file <absolute path of input file> --submit
```
withouth --submit, the jobs wont be submitted. This allows you to look at the input files.


For example
```
$ python Analysis/pat_nyx.py --input-file /projects/exasky/VizAly-Foresight/inputs/nyx/darwin_nyx_galton_wflow_test.jso --submit
```
will run CBench, PAT, and a plot generator for Cinema HTML viewers. The output will be a Cinema database.
