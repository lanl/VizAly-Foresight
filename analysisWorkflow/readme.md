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
This python script will launch batch scripts to run CBench, the halo finder, and power spectra calculation.
An example configuration file is provided, `workflow_c.ini`.

To create a Cinema database, do:
```
cd workflow_c
cp ../*.py
python create_cinema.py --input-file workflow_c.csv --debug
```


## Requirements:
workflow_b.py needs the aspw (Another Sqlite Python Wrapper). To run it on Darwin, you will need Python 3 and aspw

``` Shell Script
# load Python 3
module load anaconda/Anaconda3 openmpi/2.1.3-gcc_6.4.0 cmake/3.12.1
conda create --yes --name cbench python=3.6
source activate cbench

# install Python packages
conda install --yes numpy==1.15.4 matplotlib==3.0.2
python -m pip install apsw==3.9.2.post1
```

source the file HACC.darwin_setup to load the correct modules on Darwin.
