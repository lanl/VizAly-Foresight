#! /bin/bash

scp ${1}@${2}:${3}* .
mkdir -p output/bigcrunch_benchmark
mpirun -n 2 --allow-run-as-root ../build/bin/bigcrunch_benchmark `basename ${3}` output/bigcrunch_benchmark/result.csv 3 1 -2 0
