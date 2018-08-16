# VizAly-CBench: A Compression benchmark suite for Visualization and Analsis of Simulation Data

## Project Scope
VizAly is a general framework for **A**na**ly**sis and **Vi**suali**z**ation of simulation data. As supercomputing resources increase, cosmological scientists are able to run more detailed and larger simulations generating massive amounts of data. Analyzing these simulations with an available open-source toolkit is important for collaborative Department of Energy scientific discovery across labs, universities, and other partners. Developed software as a part of this collection include: comparing data with other existing simulations, verifying and validating results with observation databases, new halo finder algorithms, and using analytical tools to get insights into the physics of the cosmological universe. The goal of this software project is to provide a set of open-source libraries, tools, and packages for large-scale cosmology that allows scientists to visualize, analyze, and compare large-scale simulation and observational data sets. Developed software will provide a variety of methods for processing, visualization, and analysis of astronomical observation and cosmological simulation data. These tools are intended for deployment on multiple scientific computing platforms, including but not limited to personal computers, cloud computing, experimental sites (telescopes) and high-performance supercomputers.


# Installing VizAly-CBench
## Requirements:
* CMake 3.6.2 or higher
* GCC 4.9 or higher
* OpenMPI 2 or higher

## Build:
The folder scrips contains modules to load on the Cooley supercomputer at ANL and Darwin supercomputer at LANL
```
$ git clone https://github.com/lanl/VizAly-CBench.git
$ cd VizAly-CBench
$ source build.sh
```

## Running:
```
$ mpirun -np 2 ./CBench ../inputs/jesus_blosc.json
$ cat metrics_0.log
```

# Usage
CBench takes as input a json file that specifies the input parameters. The list of parameters to specity are:
* filetype (HACC or NYX)
* filename
* scalars to analyze
* log-filename
* metric-filename
* compressors to use
* metrics
Examples of input json files are in inputs folder. 

### Currently Supported file formats:
* GenericIO

### Currently Supported compressors:
* Lossless BLOSC


# Copyright and license
LANS has asserted copyright on the software package C17078, entitled Framework for Analysis and Visualization of Simulation Data.   

# Contact
Pascal Grosset, pascalgrosset@lanl.gov

[![Build Status](https://travis-ci.org/lanl/VizAly-CBench.svg?branch=master)](https://travis-ci.org/lanl/VizAly-CBench)
