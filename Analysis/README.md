# PAT: Python Analysis Toolkit

This folder contains PAT, the Python package that is used to run the Foresight workflows. The Analysis toolkit will generate a Cinema database from the input file specified.

## Prerequisites:

* Python 2.7

There is a `requirements.txt` file in the top-level of this repository. All of these packages can be installed as: `pip install -r requirements.txt`

## Running Foresight
This will run the full analysis workflow and generate a cinema database. For example
```
$ python Analysis/pat_nyx.py --input-file inputs/nyx/nyx_darwin_test.json --submit
```
will run CBench, PAT, and a plot generator for cinema. The output will be a cinema database
