# CBench
This is the compression benchmark part of Foresight. It can be run on its own or part of the Foresight workflow.

## Prerequisites:
* CMake 3.6.2 or higher
* GCC 4.9 or higher (C++11 minimum)
* OpenMPI 2 or higher


## Code Fomat
* Indentation: indentation is tabs of size 4
* Default bracket style: [Allman](https://en.wikipedia.org/wiki/Indentation_style#Allman_style)
* Varible naming: CamelCase + use meaningful variable names


## Building:
The folder **[_scripts_](scripts)** contains scripts to load modules on Cooley @ ANL and Darwin @ LANL . These build a base version of CBench:
```
$ git clone https://github.com/lanl/VizAly-Foresight.git

$ cd VizAly-Foresight
$ source buildDependencies.sh
$ source build.sh
```

## Running:
```
$ mpirun -np 2 ./CBench ../inputs/hacc/hacc_cbench_test.json
$ cat metrics_hacc_test_.csv
```
**Note:**  The above will only run a toy dataset meant for testing if the code runs. The results should **NOT** be used as an indicator for the performance of the compressors!


# Usage
CBench takes as input a json file (examples of input json files are in the **_[inputs](inputs)_** folder) that specifies the input parameters. The list of parameters to specify are:
* Filetype (HACC or NYX)
* Filename
* Scalars to analyze
* Compressor parameters
* Filename (For optional output)
* Log-filename
* Metric-filename
* Compressors to use
* Metrics

### Currently Supported file formats:
* GenericIO
* NYX(HDF5 version hdf5-1_10_3) 
* Binary

### Currently Supported compressors:
* [Lossless BLOSC](http://blosc.org/) - version 1.10.2
* [Lossy BigCrunch](https://github.com/lanl/VizAly-BigCrunch) - version 1.1
* [Lossy LossyWave](https://github.com/lanl/VizAly-LossyWave) - version 0.1
* [SZ](https://collab.cels.anl.gov/display/ESR/SZ) - version 2.0.2.2
* [zfp](https://computation.llnl.gov/projects/floating-point-compression) - version 0.5.4
* [fpzip](https://computation.llnl.gov/projects/floating-point-compression) - version 1.2.0
* [ISABELA](http://freescience.org/cs/ISABELA/ISABELA.html) - version 0.2.1

### Currently Supported metrics:
* Absolute Error
* Relative Error
* Mean Square Error (MSE)
* Peak Signal-to-Noise Ratio (PSNR)
* Memory Usage
* Compute times
