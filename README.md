# VizAly-CBench: A Compression benchmark suite for Visualization and Analsis of Simulation Data

## Project Scope
VizAly is a general framework for **A**na**ly**sis and **Vi**suali**z**ation of simulation data. As supercomputing resources increase, cosmological scientists are able to run more detailed and larger simulations generating massive amounts of data. Analyzing these simulations with an available open-source toolkit is important for collaborative Department of Energy scientific discovery across labs, universities, and other partners. Developed software as a part of this collection include: comparing data with other existing simulations, verifying and validating results with observation databases, new halo finder algorithms, and using analytical tools to get insights into the physics of the cosmological universe. The goal of this software project is to provide a set of open-source libraries, tools, and packages for large-scale cosmology that allows scientists to visualize, analyze, and compare large-scale simulation and observational data sets. Developed software will provide a variety of methods for processing, visualization, and analysis of astronomical observation and cosmological simulation data. These tools are intended for deployment on multiple scientific computing platforms, including but not limited to personal computers, cloud computing, experimental sites (telescopes) and high-performance supercomputers.



# Example

## To build:
git clone https://github.com/lanl/VizAly-CBench.git
cd VizAly-CBench
source build.sh

## Running:
mpirun -np 2 ./CBench ../inputs/jesus_blosc.json

## Verifying:
cat ../testing/data/metrics_0.log


# Copyright and license
LANS has asserted copyright on the software package C17078, entitled Framework for Analysis and Visualization of Simulation Data.   

# Contact
Pascal Grosset, pascalgrosset@lanl.gov
