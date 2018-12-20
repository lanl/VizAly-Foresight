#!/bin/bash

# Get parameter 
if [ "$#" -ne 2 ]; then 
    echo "The following paramenters are needed:"
    echo " - Path to HACC exe"
    echo " - Path to Halo input script"
    echo " - path to run path"

    exit 1
fi 


path_to_HACC_exe=$1
path_to_HALO_script=$2
run_path=$3

echo $path_to_HACC_exe
echo $path_to_HALO_script
echo $run_path


# Run Halo finder

