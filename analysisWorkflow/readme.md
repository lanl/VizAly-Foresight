The Analysis workflow will run the halo finder on specified files and then check the results of the run.

## workflow_a.py
This python script will launch batch scripts to run halo-finder on HACC raw output files. Documentation on how to build HACC and run the halo finder are located in the file HACC_readme.md
To run: python worrkflow_a.py workflow_input.json

## workflow_b.py
This python script will use the SQLite GenericIO interface to query the output files
To run: python worrkflow_a.py workflow_input.json

## Sample input file
An exmaple for the input file is workflow_input.json. workflow_input.json has two sections, A and B, marked with "Section-A" and "Section-B" respectively.

## workflow_c.py
This python script will launch batch scripts to run CBench and the halo finder.
And example configuration file is provided, `workflow_c.ini`.


## Requirements:
workflow_b.py needs the aspw (Anotger Sqlite Python Wrapper). To run it on Darwin, you will need Python 3 and aspw

``` Shell Script
#Load python 3
module load python/3.5.1

#Install apsw to user space
python -m pip install --user apsw
```

source the file HACC.darwin_setup to load the correct modules on Darwin.
