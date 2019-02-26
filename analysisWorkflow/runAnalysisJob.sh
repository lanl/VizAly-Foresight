#!/bin/bash
#SBATCH -N 1

# Go up one level before
pushd ..
projectPath=$(pwd)
popd


# load modules
source $projectPath/analysisWorkflow/HACC.darwin_setup

# go to folder 
cd $projectPath/analysisWorkflow

# Run:
python $projectPath/analysisWorkflow/raw_output_test.py
