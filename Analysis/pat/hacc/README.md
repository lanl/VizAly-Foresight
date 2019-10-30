# HACC workflow

To submit a workflow to a Slurm cluster, do:
``` Shell Script
python pat_hacc.py --input-file ../input/hacc/hacc_wflow.json --submit


python3 -m pat.hacc.workflow --input-file ../inputs/nyx/hacc_wflow.json
```


## Requirements

Here is an example of install the dependencies on Darwin:
``` Shell Script
# load Python 3
module load anaconda/Anaconda3 openmpi/2.1.3-gcc_6.4.0 cmake/3.12.1
conda create --yes --name cbench python=3.6
source activate cbench

# install Python packages
conda install --yes numpy==1.15.4 matplotlib==3.0.2
python -m pip install apsw==3.9.2.post1
```
