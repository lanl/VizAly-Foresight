#! /bin/bash

ls
pwd

source /src/env.sh

cd /src/VizAly-CBench
mkdir build
cd build
cmake ../src

make -j

./CBench ../inputs/blosc.json

cat ../testing/data/metrics.txt.log
