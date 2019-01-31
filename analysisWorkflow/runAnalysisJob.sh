#!/bin/bash
#SBATCH -N 1

projectPath=/projects/exasky/pascal-projects/VizAly-CBench

# load modules
source $projectPath/analysisWorkflow/HACC.darwin_setup

# go to folder 
cd $projectPath/analysisWorkflow

# Run:
python $projectPath/analysisWorkflow/raw_output_test.py
