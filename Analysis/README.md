# PAT: Python Analysis Toolkit

This folder contains PAT, the Python package that is used to run the Foresight workflows. The Analysis toolkit will generate a Cinema database from the input file specified.

## Prerequisites:

Minimum Requirements
* Python 3.6 or higher with (+matplotlib=3.0.2, +apsw=3.9.2, +numpy=1.15.4)
* SLURM (for job launching)

## Running Foresight
This will run the full analysis workflow and generate a cinema database. 

For example, to generate the anlysis for nyx, use the command as follows:
```
python3 <analysis_name> --input-file <absolute path of input file>
```
This will generate the scripts required to launch the full analysis.

To launch the analysis (submit the jobs to SLURM), append the '--submit' command:
```
python3 <analysis_name> --input-file <absolute path of input file> --submit
```

For example
```
$ python3 Analysis/pat_nyx.py --input-file /projects/exasky/VizAly-Foresight/inputs/nyx/darwin_nyx_galton_wflow_test.json --submit
```
will run CBench, PAT, and a plot generator for Cinema HTML viewers. The output will be a Cinema database.
