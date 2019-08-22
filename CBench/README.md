# CBench
This is the compression benchmark arm of Foresight. It can be run on its own or part of the Foresight workflow.

## Prerequisites:
* CMake 3.8.1 or higher
* GCC 6.4 or higher (C++14 minimum)
* OpenMPI 2 or higher


## Code Fomat
* Indentation: indentation is tabs of size 4
* Default bracket style: [Allman](https://en.wikipedia.org/wiki/Indentation_style#Allman_style)
* Varible naming: CamelCase + use meaningful variable names


## Building:
The folder **[_scripts_](scripts)** contains scripts to load modules on Cooley @ ANL and Darwin @ LANL . These build a base version of CBench:
```
git clone https://github.com/lanl/VizAly-Foresight.git

cd VizAly-Foresight
source buildDependencies.sh 
source build.sh
```

For more options, run buildDependencies.sh and build.sh with  "-h" 


## Running:
```
cd build/
mpirun -np 2 ./CBench ../inputs/hacc/hacc_cbench_test.json
cat metrics_hacc_test_.csv
```
**Note:**  The above will only run a toy dataset meant for testing if the code runs. The results should **NOT** be used as an indicator for the performance of the compressors!


# Usage
CBench takes as input a json file (examples of input json files are in the **_[inputs](inputs)_** folder) that specifies the input parameters.

### Currently Supported file formats:
* [GenericIO](https://trac.alcf.anl.gov/projects/genericio)
* NYX(HDF5 version hdf5-1_10_3) 
* Binary

### Currently Supported compressors:
* [Lossless BLOSC](http://blosc.org/) - version 1.10.2
* [SZ](https://collab.cels.anl.gov/display/ESR/SZ) - version 2.1.4.2
* [zfp](https://computation.llnl.gov/projects/floating-point-compression) - version 0.5.5
* [fpzip](https://computation.llnl.gov/projects/floating-point-compression) - version 1.2.0
* [ISABELA](http://freescience.org/cs/ISABELA/ISABELA.html) - version 0.2.1

### Currently Supported metrics:
* Absolute Error
* Relative Error
* Mean Square Error (MSE)
* Peak Signal-to-Noise Ratio (PSNR)
* Data Ranges (data distribution)
* Memory Usage
* Compute times
* Histogram of values
