# VizAly-Foresight: A Compression benchmark suite for Visualization and Analysis of Simulation Data

## Project Scope
VizAly is a general framework for **A**na**ly**sis and **Vi**suali**z**ation of simulation data. As supercomputing resources increase, cosmological scientists are able to run more detailed and larger simulations generating massive amounts of data. Analyzing these simulations with an available open-source toolkit is important for collaborative Department of Energy scientific discovery across labs, universities, and other partners. Developed software as a part of this collection include: comparing data with other existing simulations, verifying and validating results with observation databases, new halo finder algorithms, and using analytical tools to get insights into the physics of the cosmological universe. The goal of this software project is to provide a set of open-source libraries, tools, and packages for large-scale cosmology that allows scientists to visualize, analyze, and compare large-scale simulation and observational data sets. Developed software will provide a variety of methods for processing, visualization, and analysis of astronomical observation and cosmological simulation data. These tools are intended for deployment on multiple scientific computing platforms, including but not limited to personal computers, cloud computing, experimental sites (telescopes) and high-performance supercomputers.

Foresight has three components:
* CBench: the compression benchmark arm of Foresight desighned to run at scale on supercomputers
* PAT: Python Analysis Toolkit, which contains a bunch of untilities to speed up analysis and plotting of results
* Cinema: The tool used to visualize the results of this project

## Current Visualization results
* Link to the live cinema database with results : [https://lanl.github.io/VizAly-Foresight/](https://lanl.github.io/VizAly-CBench/)
* Cinema comparison is at: https://lanl.github.io/VizAly-CBench/cinema_compare/

# Building VizAly-Foresight
The default master branch should always point to the latest working version. However, for more stable releases, you should checkout the latest tag release.

## Prerequisites:
* CMake 3.8.1 or higher
* GCC 6.4 or higher (C++11 minimum)
* OpenMPI 2 or higher

## Building:
The folder **[_scripts_](scripts)** contains scripts to load modules on Cooley @ ANL and Darwin @ LANL . These build a base version of Foresight:
```
$ git clone https://github.com/lanl/VizAly-Foresight.git

$ cd VizAly-Foresight
$ source buildDependencies.sh
$ source build.sh
```

To build a master (all) version of Foresight, run the following scripts:
```
$ git clone https://github.com/lanl/VizAly-Foresight.git

$ cd VizAly-Foresight
$ source buildAllDependencies.sh
$ source buildAll.sh
```

## Running:
```
$ mpirun -np 2 ./CBench ../inputs/hacc/HACC_all.json
$ cat metrics_HACC_all_
```

**Note:**  The above will only run a toy dataset meant for testing if the code runs. The results should **NOT** be used as an indicator for the performance of the compressors!

### Tools:
[tools/plotting](tools/plotting) contains a python program for easy graphing on the csv results

# Usage
Foresight takes as input a json file (examples of input json files are in the **_[inputs](inputs)_** folder) that specifies the input parameters. The list of parameters to specify are:
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

# Development
For information on how to add new compressors and/or metrics, look at the [readme in CBench/compressors](CBench/compressors/readme.md) and [CBench/metrics](CBench/metrics) respectively.

For information on Travis CI and Docker image, look at the [travis folder](testing/travis) 

# Developers
* Chris Biwer
* Pascal Grosset
* Jesus Pulido
* Hoby Rakotoarivelo

[![Build Status](https://travis-ci.org/lanl/VizAly-Foresight.svg?branch=pat_refactor_pascal;)](https://travis-ci.org/lanl/VizAly-Foresight)


# Copyright and License
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Triad National Security, LLC.
All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
