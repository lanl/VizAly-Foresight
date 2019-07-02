# PAT: This folder contains PAT, the Python Analysis Toolkit that is used to run the foresight workflow. The Analysis toolkit will generate a cinema database from the input file specified

## Prerequisites:
* Python 2
* numpy
* pandas
* matplotlib

All of these packages can be installed as python pip install <package name>

## Running Foresight
This will run the full analysis workflow and generate a cinema database. For example
```
$ python Analysis/pat_nyx.py --input-files inputs/nyx/nyx_darwin_test.json --submit
```
will run CBench, PAT, and a plot generator for cinema. The output will be a cinema database
