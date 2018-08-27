# VizAly-CBench: A Compression benchmark suite for Visualization and Analysis of Simulation Data

## Project Scope
VizAly is a general framework for **A**na**ly**sis and **Vi**suali**z**ation of simulation data. As supercomputing resources increase, cosmological scientists are able to run more detailed and larger simulations generating massive amounts of data. Analyzing these simulations with an available open-source toolkit is important for collaborative Department of Energy scientific discovery across labs, universities, and other partners. Developed software as a part of this collection include: comparing data with other existing simulations, verifying and validating results with observation databases, new halo finder algorithms, and using analytical tools to get insights into the physics of the cosmological universe. The goal of this software project is to provide a set of open-source libraries, tools, and packages for large-scale cosmology that allows scientists to visualize, analyze, and compare large-scale simulation and observational data sets. Developed software will provide a variety of methods for processing, visualization, and analysis of astronomical observation and cosmological simulation data. These tools are intended for deployment on multiple scientific computing platforms, including but not limited to personal computers, cloud computing, experimental sites (telescopes) and high-performance supercomputers.


# Building VizAly-CBench
The default master branch should always point to the latest working version. However, for more stable releases, you should check out the latest tag release. The current latest is v1.0

## Prerequisites:
* CMake 3.6.2 or higher
* GCC 4.9 or higher (C++11 minimum)
* OpenMPI 2 or higher

## Building:
The folder **_scripts_** contains scripts to load modules on Cooley @ ANL and Darwin @ LANL
```
$ git clone https://github.com/lanl/VizAly-CBench.git

$ cd VizAly-CBench
$ source build.sh
```

## Running:
```
$ mpirun -np 2 ./CBench ../inputs/all.json
$ cat metrics
```

# Usage
CBench takes as input a json file (examples of input json files are in the **_inputs_** folder) that specifies the input parameters. The list of parameters to specify are:
* Filetype (HACC or NYX)
* Filename
* Scalars to analyze
* Compressor parameters
* Log-filename
* Metric-filename
* Compressors to use
* Metrics

### Currently Supported file formats:
* GenericIO

### Currently Supported compressors:
* [Lossless BLOSC](http://blosc.org/)
* [Lossy BigCrunch](https://github.com/lanl/VizAly-BigCrunch)
* [SZ](https://collab.cels.anl.gov/display/ESR/SZ)

### Currently Supported metrics:
* Absolute Error
* Relative Error
* MSE

# Development
For information on how to add new compressors and/or metrics, look at the [tutorial folder](https://github.com/lanl/VizAly-CBench/tree/master/tutorial)

For information on how to develop the Travis CI and Docker image, look at the [travis folder](https://github.com/lanl/VizAly-CBench/tree/master/testing/travis)

# Copyright and license
LANS has asserted copyright on the software package C17078, entitled Framework for Analysis and Visualization of Simulation Data.   

# Developers
* Chris Biwer
* Pascal Grosset
* Jesus Pulido

[![Build Status](https://travis-ci.org/lanl/VizAly-CBench.svg?branch=master)](https://travis-ci.org/lanl/VizAly-CBench)
