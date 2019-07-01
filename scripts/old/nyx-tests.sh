#! /bin/bash

num_ranks=4
project="$(pwd)/.."
build_dir=$project/build
input_dir=$project/inputs
script_dir=$project/scripts

clean() {
  printf "Cleaning build directory ... "
  cd "$build_dir" && \
  rm -f *.h5 *.log core.* *.hdf5 metrics_nyx_* && \
  printf "done\n"
}

build() {
  cd "$project" && source buildAll.sh
}

run() {
  echo "Running testcases ..."
  cd "$build_dir"

  for testcase in $input_dir/tests/*.json; do
    if [ -e "$testcase" ]; then
      echo "mpirun -np $num_ranks CBench $testcase" 
      mpirun -np "$num_ranks" CBench "$testcase" 
    fi
  done
}

#source clean.sh && cd .. && source build.sh && mpirun -np 4 CBench ../inputs/nyx_attr_test.json 
clean && build && run

