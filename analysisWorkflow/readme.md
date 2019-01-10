The Analysis workflow will run the halo finder on specified files and then check the results of the run.

## workflow_a.py
This python script will launch batch scripts to run halo-finder on HACC raw output files

## workflow_b.py
This python script will use the SQLite GenericIO interface to query the output files

All the input are in the file workflow_input.json


To run this on Darwin, you will need Python 3 and aspw

``` Shell Script
#Load python 3
module load python/3.5.1

#Install apsw to user space
python -m pip install --user apsw
```
