#!/bin/bash

# Run compressors benchmarks on NYX dataset
# Author: Hoby Rakotoarivelo

numNodes=16
ranksPerNode=8
project="$EXASKY_HOME/cbench"
#compressors=( bigcrunch fzip isabela sz zfp )
compressors=( all )

setup() {
  module purge
  module load python/2.7.3
  module load gcc/6.4.0
  module load openmpi/2.1.3-gcc_6.4.0
  module load cmake
 
  export LD_LIBRARY_PATH=/projects/exasky/fftw-3.3.8/install/lib:$LD_LIBRARY_PATH
  export FFTW_MAJOR_VERSION=3
  
  export NYX_PLATFORM="Darwin"
  export NYX_OBJDIR="${NYX_PLATFORM}"
  
  export NYX_CC="gcc"
  export NYX_CXX="g++"
  export NYX_CFLAGS="-O3 -g -fopenmp -std=c++11"
  export NYX_CXXFLAGS="-O3 -g -fopenmp -std=c++11"
  export NYX_LDFLAGS="-lm -fopenmp"
  
  export NYX_MPI_CC="mpicc"
  export NYX_MPI_CXX="mpicxx"
  export NYX_MPI_CFLAGS="-O3 -std=gnu99 -g -fopenmp -std=c++11"
  export NYX_MPI_CXXFLAGS="-O3 -g -Wno-deprecated -fopenmp -fPIC -std=c++11"
  export NYX_MPI_LDFLAGS="-lm -fopenmp"
  
  export GIO_MPICXX="${NYX_MPI_CXX}"
}

submit() {

  for compressor in "${compressors[@]}"; do
    #input="../inputs/nyx_therma_l_512_${compressor}.json"
    input="../inputs/nyx/nyx_thermal_reduced.json"
    script="nyx-${compressor}.sh"
    #config="$EXASKY_ROOT/VizAly-CBench/scripts/VizAly-CBench.bash.darwin"
    binary="$project/build/CBench"
    
    echo "#!/bin/bash"                               > $script
    echo "#SBATCH -N $numNodes"                     >> $script
    echo "#SBATCH --ntasks-per-node $ranksPerNode"  >> $script
    echo "#SBATCH -p scaling"                       >> $script
    echo "source $config"                           >> $script
    echo "cd $project/build"                        >> $script
    echo "mpirun $binary $input"                    >> $script
    
    sbatch $script
  done
}

clean() {
  # remove temp scripts
  for compressor in "${compressors[@]}"; do
    script="nyx-${compressor}.sh"
    [ -f $script ] && rm $script
  done	

  # sort output files
  #source "$project/build/clean.sh"	
}

# main
setup && submit && clean
  
