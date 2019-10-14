# VizAly-Foresight: A Compression benchmark suite for Visualization and Analysis of Simulation Data

## Project Scope
VizAly is a general framework for **A**na**ly**sis and **Vi**suali**z**ation of simulation data. As supercomputing resources increase, cosmological scientists are able to run more detailed and larger simulations generating massive amounts of data. Analyzing these simulations with an available open-source toolkit is important for collaborative Department of Energy scientific discovery across labs, universities, and other partners. Developed software as a part of this collection include: comparing data with other existing simulations, verifying and validating results with observation databases, new halo finder algorithms, and using analytical tools to get insights into the physics of the cosmological universe. The goal of this software project is to provide a set of open-source libraries, tools, and packages for large-scale cosmology that allows scientists to visualize, analyze, and compare large-scale simulation and observational data sets. Developed software will provide a variety of methods for processing, visualization, and analysis of astronomical observation and cosmological simulation data. These tools are intended for deployment on multiple scientific computing platforms, including but not limited to personal computers, cloud computing, experimental sites (telescopes) and high-performance supercomputers.

Foresight has three components:
* CBench: the compression benchmark arm of Foresight desighned to run at scale on supercomputers
* PAT: Python Analysis Toolkit, which contains a bunch of untilities to speed up analysis and plotting of results
* Visualization: [Cinema](https://cinemascience.github.io/) is used to visualize the results of this project 

## Current Visualization results
* Link to the live cinema database with results : [https://lanl.github.io/VizAly-Foresight/](https://lanl.github.io/VizAly-Foresight/)
* Cinema comparison is at: https://lanl.github.io/VizAly-Foresight/cinema_compare/

# Building VizAly-Foresight
The default master branch should always point to the latest working version. However, for more stable releases, you should checkout the latest tag release.

## Prerequisites:
CBench
* CMake 3.12 or higher
* GCC 6.4 or higher (C++14 minimum)
* OpenMPI 2 or higher

PAT
* Python 3.6 or higher with (+matplotlib=3.0.2, +apsw=3.9.2, +numpy=1.15.4)
* SLURM (for job launching)

Visualization
* Mozilla Firefox Recommended for Cinema Explorer

## Building:
The folder **[_scripts_](scripts)** contains scripts to load modules on Cooley @ ANL, Cori @ NERSC, and Darwin @ LANL . These build a base version of Foresight:
```
git clone https://github.com/lanl/VizAly-Foresight.git

cd VizAly-Foresight
source scripts/<Name of the environment> # that sets up the environment
source buildDependencies.sh
source build.sh
```

## Running Foresight
This will set up and run the full analysis workflow on SLURM and generate a cinema database, the command is as follows:
```
cd Analysis
python3 -m <name_of_analysis> --input-file <path to input JSON file>
```
For example, to run the NYX analysis, the command is:
```
python3 -m pat.nyx.workflow --input-files inputs/nyx/nyx_darwin_test.json
```

The outputs/logs will be generated in a folder using "project-home" and "wflow-path" in the input JSON file. The folder 
cinemaDB.cdb will contain the graphed outputs.


## Running CBench as stand-alone:
```
mpirun -np 2 build/CBench ../inputs/hacc/hacc_cbench_test.json
cat metrics_hacc_test_.csv
```
Alternatively, CBench can also be run throughj foresight as follows:
```
cd Analysis
python3 -m <name_of_analysis> --input-file <path to input JSON file> --cbench
```

**Note:**  The above will only run a toy dataset meant for testing if the code runs. The results should **NOT** be used as an indicator for the performance of the compressors!


### Provenance
The folder in cinemaDB.cdb will also contain a wflow.json that will point to a git-hash tag that the code was ran on. To replicate the run, git checkout <git-tag> to get the exact same code.


# Usage
Foresight takes as input a json file (examples of input json files are in the **_[inputs](inputs)_** folder) that specifies the input parameters. 

### Currently Supported file formats:
* GenericIO
* NYX(HDF5 version hdf5-1_10_3) 
* Binary

### Currently Supported compressors:
* [Lossless BLOSC](http://blosc.org/) - version 1.10.2
* [SZ](https://collab.cels.anl.gov/display/ESR/SZ) - version 2.1.4.2
* [zfp](https://computation.llnl.gov/projects/floating-point-compression) - version 0.5.5
* [fpzip](https://computation.llnl.gov/projects/floating-point-compression) - version 1.2.0
* [ISABELA](http://freescience.org/cs/ISABELA/ISABELA.html) - version 0.2.1

The following experimental **GPU** compressors are also supported:
* SZ - custom version
* [zfp](https://github.com/LLNL/zfp.git) - version 0.5.4

### Currently Supported metrics:
* Absolute Error
* Relative Error
* Mean Square Error (MSE)
* Peak Signal-to-Noise Ratio (PSNR)
* Data Ranges (Data distribution)
* Memory Usage
* Compute times

# Development
For information on how to add new compressors and/or metrics, look at the [readme in CBench/compressors](CBench/compressors/readme.md) and [CBench/metrics](CBench/metrics) respectively. To add new analysis routines, look at the [readme in Analysis/ folder](Analysis/README.md).

For information on Travis CI and Docker image, look at the [travis folder](testing/travis) 

# Developers
* Chris Biwer
* Pascal Grosset
* Sian Jin
* Jesus Pulido
* Hoby Rakotoarivelo

[![Build Status](https://travis-ci.org/lanl/VizAly-Foresight.svg?branch=master;)](https://travis-ci.org/lanl/VizAly-Foresight)


# Copyright and License
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Triad National Security, LLC.
All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
