#! /bin/bash

set -e

# get GCC compiler
source /src/env.sh

# build
cd /src/VizAly-CBench
projectPath=$(pwd)

source buildDependencies.sh	# build dependencies
source build.sh				# build the code

# run example
mpirun -np 2 --allow-run-as-root ./CBench ../inputs/HACC_all.json

# view output
cat metrics_HACC_all_
