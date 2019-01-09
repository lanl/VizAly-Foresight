# HACC
This document contains information to build and run HACC on the Darwin supercomputer.

## i) Getting HACC
If you have access to the HACC repo at Argonne, checkout the code to get the latest version:
~~~~
svn checkout https://svn.alcf.anl.gov/repos/DarkUniverse/hacc/trunk
~~~~

That should create a folder called trunk. 

## ii) Building HACC
We first need to load a suitable environment for HACC:
~~~~
source /projects/exasky/HACC.darwin_setup 
~~~~


where HACC.darwin_setup is:
~~~~
module purge
module load python/2.7.3
module load gcc/6.4.0
module load openmpi/2.1.3-gcc_6.4.0
module load cmake

export FFTW_MAJOR_VERSION=3

export HACC_PLATFORM="Darwin"
export HACC_OBJDIR="${HACC_PLATFORM}"

export HACC_CFLAGS="-O3 -g -fopenmp -std=c++11"
export HACC_CC="gcc"

export HACC_CXXFLAGS="-O3 -g -fopenmp -std=c++11"
export HACC_CXX="g++"

export HACC_LDFLAGS="-lm -fopenmp"

export HACC_MPI_CFLAGS="-O3 -std=gnu99 -g -fopenmp -std=c++11"
export HACC_MPI_CC="mpicc"


export HACC_MPI_CXXFLAGS="-O3 -g -Wno-deprecated -fopenmp -fPIC -std=c++11"
export HACC_MPI_CXX="mpicxx"

export HACC_MPI_LDFLAGS="-lm -fopenmp"

## Propagate GenericIO env variables
export GIO_MPICXX="${HACC_MPI_CXX}"
~~~~

After loading a suitable environment, we can start building HACC:
~~~~
cd trunk
make -j16
~~~~

The compilation output will be in the “Darwin”subfolder and subfolders mpi and frontend under which there will be lib, bin, and other folders.

## iii) Running HACC
The following are required to run the halo finder:
* Executable hacc_pm
* Input configuration parameters indat.params
* Analysis parameter analysisdat
* Another input file parameter file cmbM000.tf

The following setup is needed to run HACC; and output folder and some input parameters
~~~~
cd $PROJECT
mkdir run
mkdir output
mkdir inputs
cd run 
~~~~

where $PROJECT in this case is /projects/exasky/HACC/. Use the file indat_one.params to do a test run as follows:
~~~~
cd $PROJECT/trunk
mpirun -np 8 $PROJECT/trunk/Darwin/mpi/bin/hacc_pm -n inputs/indat_one.params
~~~~

This should produce a small test file. To run a test on multiple nodes, a batch job is needed. An example script (for 256 ranks on 16 nodes with 16 mpi ranks per node) is:
Jobscript.sh
~~~~
#!/bin/bash
#SBATCH -N 16
#SBATCH --ntasks-per-node 16
#SBATCH -p scaling


# load modules
source /projects/exasky/HACC.darwin_setup 


# go to folder 
cd /projects/exasky/HACC/run/


# Run:
mpirun /projects/exasky/HACC/trunk/Darwin/mpi/bin/hacc_pm -n inputs/indat.params
~~~~

Note:
* If there are errors associated to loopback, add “export OMPI_MCA_btl_tcp_if_exclude=virbr0,lo” to the jobscript.

Launch as follows:
~~~~
sbatch jobscript.sh
~~~~


Some example input files are here /projects/exasky/HACC/run/inputs on Darwin and available here.
analysisdat
~~~~
0.168   # linking length
1.0e5   # minimum mass (M_solar/h)
40      # minimum number of particles
5.0     # overload length (Mpc/h)
22500   # Nskewers
4096    # Npixels
0.325   # h (Mpc/h): pixel smoothing length
32      # Nmesh<  L/h
1.0     # dist. conversion factor
1       # output FOF halo property summary
1       # enable outputs with particles (may be large; controlled by options below)
1       # output SOD halo property summary
0       # output subhalo property summary
1.0     # rho_c conversion factor
1.0     # SOD-mass conversion factor
1000    # the minimum number of particles for which to compute SOD properties
1       # run the MBP algorithm for FOF centers
0       # run the MCP algorithm for FOF centers
0       # use the minimum potential array
64      # number of neighbors for calculating density
20      # number of close neighbors used in subgrouping
20      # minimum particles in a subhalo
20000   # smallest FOF halo to have subfinding run on
1.0     # factor for cut/grow criteria
0.0     # factor for Poisson noise significance
1       # min. number of neighbors for linking
10000  # output particles in halos with particles larger than this number of particles
0.01    # fraction of all particles in halos meeting the halo size cut to output
5       # the minimum number of particles per halo in the fractional output
1       # output all particle tags
0       # output all particles
~~~~

indat.params 
~~~~
################################################################################
# Header version information
################################################################################
HACC_HEADER_VERSION 1.0.0

################################################################################
# Cosmological Parameters
# Length scales are measured in Mpc/h
# OMEGA_CDM and OMEGA_NU given for Omega_cdm and Omega_nu (no $h^2$)
# DEUT=Omegab*h^2 
# HUBBLE: Hubble constant/100 km/s/Mpc
# SS8: target value for sigma_8
# NS: index of the primordial power spectrum
# W_DE: constant dark energy equation of state
# Currently flat Universe only
################################################################################
OMEGA_CDM 0.220
DEUT 0.02258
OMEGA_NU 0.0
HUBBLE 0.71
SS8 0.8
NS 0.963
W_DE -1.0
WA_DE 0.0

################################################################################
# Initializer Set-up and read-ins
# ZIN: Starting redshift
# USE_WHITE_NOISE_INIT: YES: real space, NO: k space
# input type: INIT|RECORD|BLOCK|COSMO|RESTART
# INIT: generates internal initial conditions, rest if for read-ins
# distrib. type: ROUND_ROBIN|ALL_TO_ALL|ONE_TO_ONE|restart_step
#                (ignored if INPUT_TYPE is INIT)
# ROUND_ROBIN indicates particles must be looked at by all processors
# ONE_TO_ONE indicates that particles physically reside on matching processor
# ALL_TO_ALL improved ROUND_ROBIN
# For restart: specify time step and modify INPUT_BASE_NAME
# TRANS: Transfer function: Read in CAMB file (specify name in INPUT_BASE_NAME) 
#        or internal TF (KH, HS, PD, BBKS)
################################################################################
Z_IN 50.0
USE_WHITE_NOISE_INIT YES
TRANS CMB
INPUT_BASE_NAME inputs/cmbM000.tf
INPUT_TYPE INIT
DISTRIBUTE_TYPE 170

################################################################################
# Outputs for initial conditions, alive particles, some analysis and restarts, refreshes
# WRITE_IC: write initial condition, format will be the same as for all outputs
# USE_MPI_IO: YES=one large file, NO=one file per rank in cosmo format
# REFRESH: takes either explicit time steps or file name with time steps specified,
#          same is true for all other outputs besides FINAL_GRID and VIS_SLAB
#          which only happen at the end if commented in
# SMALL_DUMP: prints all particles from rank 0
# OUTPUT_FRACTION: fraction of particles in alive dumps
# VIS_STEP: prints uniform grid of full simulation in Insley format
# FINAL_GRID_OUTPUT: ascii file! prints grid at last time step, only for small runs
# VIZ_SLAB: prints slice of final grid in Insley format
################################################################################
WRITE_IC NO
USE_MPI_IO YES
OUTPUT_BASE_NAME output/_m000
REFRESH 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400 420 440 460 480
ALIVE_DUMP 
OUTPUT_FRACTION 0.01
FULL_ALIVE_DUMP 0 50 100 150 200 250 300 350 400 450 499
SMALL_DUMP 
RESTART_DUMP 
PK_DUMP 
#VIZ_STEP 30
FINAL_GRID_OUTPUT NO
VIZ_SLAB 

################################################################################
# Code parameters I: essential parameters to specify the run/resolution
# ISEED: random number for realizatio
# NG: number of grid points (1d), NP: number of particles (1d)
# RL: physical box size [h^(-1)Mpc]
# Z_FIN: final redshift
# Timestepper: N_STEPS: number of PM steps, N_SUB: number of sub-cycles (3-5)
# OL: PM overload and refresh: 8 Mpc good choice, depends on overall volume
# RSM: Tree smoothing scale, fraction of PM grid size
# max RCB tree particles per leaf, ~100 optimal for BG/Q, ~24 for X86
################################################################################
I_SEED 7828967
NG 512
NP 512
RL 256.0
Z_FIN 0.0
N_STEPS 500
N_SUB 2
OL 8.0
RSM 0.007
RCB_TREE_PPN 24

################################################################################
# Code parameters II: specifications for tree/memory etc., mostly unchanged
# CM_SIZE: chaining mesh size, 3.12 good choice, should not be smaller
# OPENING_ANGLE: tree code
# EDGE: ???
# alpha, power of scale factor in timestepping
# TOPOLOGY: allows user to pick 3d decomposition (=ranks), if commented out, 
#           machine will pick it for you
# BIGCHUNK: memory management, percent extra to allocate in bigchunk
# USE_ALLV: important for N-to-N write, will not work on Mira at scale
################################################################################
CM_SIZE 3.12	
OPENING_ANGLE 0.1
EDGE 3.2
ALPHA 1.0
#TOPOLOGY 16x16x8
USE_BIGCHUNK YES
BIGCHUNK_EXTRA_FACTOR 5
MEMORY_PADDING_DISPLACEMENT 8.0
USE_MONOPOLE_RCB_TREE YES
USE_ALLTOALLV YES
USE_POLY YES

################################################################################
# Analysis
# Config file for analysis
################################################################################
ANALYSIS_CONFIG inputs/analysisdat
ENABLE_HALO_OUTPUTS YES
STATIC_DUMP 499
~~~~

## iv) Running HACC Halo Finder Offline
The following are required to run the halo finder:
* Executable cosmotools-driver
* input data m000.full.mpicosmo
* Timestep timesteps.dat
* halo configuration parameters cosmotools-config.dat
* Hacc configuration parameters indat.params

It is run as follows:
~~~~
$mpirun <path_to_exec>/cosmotools-driver --config <path_to_halo_inputs>/cosmotools-config.dat --timesteps <path_to_halo_inputs>/timesteps.dat --prefix <path_to_input_data>/m000.full.mpicosmo 
<path_to_halo_inputs>/indat.params
~~~~

Example run line:
~~~~
$mpirun /projects/exasky/HACC/trunk/Darwin/mpi/bin/cosmotools-driver --config /projects/exasky/HACC/run/haloInputs/cosmotools-config.dat --timesteps /projects/exasky/HACC/run/haloInputs/timesteps.dat --prefix /projects/groups/vizproject/dssdata/cosmo/Argonne_L360_HACC001/STEP499/m000.full.mpicosmo /projects/exasky/HACC/run/haloInputs/indat.params
~~~~


Example Job script:
~~~~
#!/bin/bash
#SBATCH -N 16
#SBATCH --ntasks-per-node 16
#SBATCH -p scaling

# load modules
source /projects/exasky/HACC.darwin_setup 

# go to folder 
cd /projects/exasky/HACC/run/

# Run:
mpirun /projects/exasky/HACC/trunk/Darwin/mpi/bin/cosmotools-driver --config /projects/exasky/HACC/run/haloInputs/cosmotools-config.dat --timesteps /projects/exasky/HACC/run/haloInputs/timesteps.dat --prefix /projects/groups/vizproject/dssdata/cosmo/Argonne_L360_HACC001/STEP499/m000.full.mpicosmo /projects/exasky/HACC/run/haloInputs/indat.params
~~~~

The following files are required:

timesteps.dat
~~~~
499
~~~~

indat.params 
~~~~
#####################################################################
# Header version information
#####################################################################
HACC_HEADER_VERSION 1.0.0
#####################################################################
# Cosmological Parameters
# Length scales are measured in Mpc/h
# OMEGA_CDM and OMEGA_NU given for Omega_cdm and Omega_nu (no $h^2$)
# DEUT=Omegab*h^2
# HUBBLE: Hubble constant/100 km/s/Mpc
# SS8: target value for sigma_8
# NS: index of the primordial power spectrum
# W_DE: constant dark energy equation of state
# Currently flat Universe only
#####################################################################
OMEGA_CDM 0.259
DEUT 0.022
OMEGA_NU 0.0
HUBBLE 0.676
SS8 0.817
NS 0.963
W_DE -0.92
WA_DE -0.4
T_CMB 0.0
N_EFF_MASSLESS 0.0
N_EFF_MASSIVE 0.0
#####################################################################
# Initializer Set-up and read-ins
# ZIN: Starting redshift
# USE_WHITE_NOISE_INIT: YES: real space, NO: k space
# input type: INIT|RECORD|BLOCK|COSMO|RESTART
# INIT: generates internal initial conditions, rest if for read-ins
# distrib. type: ROUND_ROBIN|ALL_TO_ALL|ONE_TO_ONE|restart_step
#            	(ignored if INPUT_TYPE is INIT)
# ROUND_ROBIN indicates particles must be looked at by all processors
# ONE_TO_ONE indicates that particles physically reside on matching processor
# ALL_TO_ALL improved ROUND_ROBIN
# For restart: specify time step and modify INPUT_BASE_NAME
# TRANS: Read in CAMB transfer function file (specify name in INPUT_BASE_NAME)
#    	or internal TF (KH, HS, PD, BBKS)
#####################################################################
Z_IN 200.0
USE_WHITE_NOISE_INIT YES
USE_CBRNG YES
TRANS CMB
INPUT_BASE_NAME cmbM001.tf
INPUT_TYPE INIT
DISTRIBUTE_TYPE LAST
MAX_MINUTES 120
###################################################################### 
# Outputs for initial conditions, alive particles, some analysis and restarts, refreshes
# WRITE_IC: write initial condition, format will be the same as for all outputs
# USE_MPI_IO: YES=one large file, NO=one file per rank in cosmo format
# REFRESH: takes either explicit time steps or file name with time steps specified,
#      	same is true for all other outputs besides FINAL_GRID and VIS_SLAB
#      	which only happen at the end if commented in
# SMALL_DUMP: prints all particles from rank 0
# OUTPUT_FRACTION: fraction of particles in alive dumps
# VIS_STEP: prints uniform grid of full simulation in Insley format
# FINAL_GRID_OUTPUT: ascii file! prints grid at last time step, only for small runs
# VIZ_SLAB: prints slice of final grid in Insley format
#####################################################################
WRITE_IC NO
USE_MPI_IO YES
OUTPUT_BASE_NAME ../../output/m001
REFRESH 10 20 
ALIVE_DUMP
OUTPUT_FRACTION 0.01
FULL_ALIVE_DUMP 97
SMALL_DUMP
RESTART_DUMP
PK_DUMP 
VIZ_STEP 
FINAL_GRID_OUTPUT NO
VIZ_SLAB
#####################################################################
# Code parameters I: essential parameters to specify the run/resolution
# ISEED: random number for realizatio
# NG: number of grid points (1d), NP: number of particles (1d)
# RL: physical box size [h^(-1)Mpc]
# Z_FIN: final redshift
# Timestepper: N_STEPS: number of PM steps, N_SUB: number of sub-cycles (3-5)
# OL: PM overload and refresh: 8 Mpc good choice, depends on overall volume
# RSM: Tree smoothing scale, fraction of PM grid size
# max RCB tree particles per leaf, ~100 optimal for BG/Q, ~24 for X86
#####################################################################
I_SEED 99847958
NG 1344
NP 1344
RL 250.0
Z_FIN 0.0
N_STEPS 500
N_SUB 1
OL 8.0
RSM 0.01
RCB_TREE_PPN 128
#####################################################################
# Code parameters II: specifications for tree/memory etc., mostly unchanged
# CM_SIZE: chaining mesh size, 3.12 good choice, should not be smaller
# OPENING_ANGLE: tree code
# EDGE: ???
# alpha, power of scale factor in timestepping
# TOPOLOGY: allows user to pick 3d decomposition (=ranks), if commented out,
#       	machine will pick it for you
# BIGCHUNK: memory management, percent extra to allocate in bigchunk
# USE_ALLV: important for N-to-N write, will not work on Mira at scale
#####################################################################
CM_SIZE 8.0
USE_CHAINING_MESH YES
CHAINING_MESH_THREADS 8
CHAINING_MESH_PER_SUBCYCLE YES
RCB_TREE_EXTRA_LEVELS 3
OPENING_ANGLE 0.1
EDGE 3.2
ALPHA 1.0
#TOPOLOGY 32x32x16
USE_BIGCHUNK YES
BIGCHUNK_EXTRA_FACTOR 20
MEMORY_PADDING_DISPLACEMENT 10.0
USE_MONOPOLE_RCB_TREE YES
USE_ALLTOALLV YES
#####################################################################
# Analysis
# Config file for analysis
#####################################################################
COSMOTOOLS ON
COSMOTOOLS_CONFIG cosmotools-config.dat
~~~~



cosmotools-config.dat 
~~~~
##==============================================================##
COSMOTOOLS CONFIGURATION
##============================================================= ##

# Set the version of the configuration file (for backwards compatibility etc.)
VERSION 1.0

# Visualization Parameters
VISUALIZATION  NO
VIZ_SERVER 	127.0.0.1
VIZ_PORT   	2222

## Frequency at which to update visualization, e.g., every 20 time-steps.
VIZ_FREQUENCY 20

## Frequency at which this configuration file will be read in to update
## any in situ analysis parameters.
CONFIG_REFRESH_FREQUENCY -1

# Enable/Disable the tools that will be used at each time-step.
# The name of the analysis tool can be anything as long as it has no spaces.
ANALYSISTOOL HALOFINDER 	YES
ANALYSISTOOL CHECKPOINT_IO  NO

##============================================================= ##
##  IN SITU ALGORITHM PARAMETERS
##============================================================= ##

##============================================================= ##
##  HALOFINDER PARAMETERS
##============================================================= ##
SECTION HALOFINDER

## Name of the internal AnalysisTool instance used
INSTANCE_NAME LANLHALOFINDER

##-----------------------------------------------------| Framework Parameters |

## Frequency type: EXPLICIT(0) or IMPLICIT(1) or REDSHIFTS(2)
FREQUENCY_TYPE 0
EXPLICIT_TIMESTEPS 499
IMPLICIT_TIMESTEPS
WRITE_OUTPUT  YES
BASE_OUTPUT_FILE_NAME  ../run/output/analysis/Halos/b0168/m001

## Indicate whether halos will be visualized 
VISIBLE YES

##-----------------------------------------------------| Algorithm Parameters |

COMPUTE_ELLIPTICITY NO  # Enable halo ellipticity output (currently will report in both FOF and SOD halo property files

#-------------------------------------| FOF Parameters |-----------------------

LINKING_LENGTH  	0.168 # The linking length used for FOF algorithm
MINIMUM_MASS    	1.0e5 # The minimum mass (M_solar/h)
FOF_PMIN        	40	# Minimum number of particles to consider as a halo
OVERLOAD_LENGTH 	6.0   # (Mpc/h)
DIST_CONVERT_FACTOR 1.0   # positions are multiplied by this factor
NEIGH_MIN       	1 	# min number of neighbors for linking (not used)
SMOOTHING_LENGTH 0.001	# smoothing length for center finder

#---------------------------| Center Finder Parameters |-----------------------

COMPUTE_FOF_CENTERS YES   # enables finding of halo centers
USE_MBP_FINDER  	YES   # run the MBP algorithm for FOF centers
USE_MCP_FINDER  	NO	# run the MCP algorithm for FOF centers
USE_HIST_FINDER 	NO	# run the histogram MCP algorithm for FOF centers
MAX_FOR_CENTER_FINDING 0 # run center finder for halos with this maximum number of particles


#---------------------------|  Parameters |-----------------------

COMPUTE_ELLIPTICITY NO  # Enable halo ellipticity output (currently will report in both FOF and SOD halo property files)


#---------------------------| Core Finder Parameters |-------------------------
NUM_CORE_SIZE 20
CORE_HALO_SIZE 100
ACCUMULATE_CORE_NAME ../run/output/analysis/Halos/b0168/m001

#-------------------------------------| SOD Parameters |-----------------------

COMPUTE_SOD_HALOS	YES	# enables SOD halos (requires FOF centers)
RHO_C_CONVERT_FACTOR 1.0   # rho_c conversion factor
SOD_MASS_CONVERSION  1.0   # SOD-mass conversion factor
SOD_PMIN         	500  # min. num particles for which SOD properties
RHO_RATIO        	200.0 # overdensity value
CDELTA_MINIMUM_PARTICLES 500 # min. num particles for concentration measurement
#---------------------------------| Subhalo Parameters |-----------------------

COMPUTE_SUBHALOS   NO	# enables/disables sub-halofinding
USE_SUBFIND    	NO	# enables/disables SUBFIND
NEIGHS_SPH_DENSITY 64	# number of neighbors for computing density
NUM_SUBHALO_NEIGHS 20	# number of close neighbors used in sub-grouping
SUBHALO_PMIN   	20	# minimum particles in a sub-halo
FOF_MIN_SIZE   	5000  # size of the smallest FOF halo to have sub-finding
ALPHA          	1.0   # factor for cut/grow criteria
BETA           	0.0   # factor for Poisson noise significance

#----------------------------------| Output Parameters |-----------------------

OUTPUT_FOF_HALO_SUMMARY YES   # generate a FOF halo property summary
OUTPUT_FOF_COM      	NO	# output center of mass in FOF halo summary
OUTPUT_SOD_HALO_SUMMARY YES	# generate a SOD halo property summary
OUTPUT_SOD_COM      	NO	# output center of mass in SOD halo summary
OUTPUT_SUBHALO_SUMMARY  NO	# generate a sub-halo property summary
OUTPUT_SUBHALO_COM  	NO	# output center of mass in sub-halo summary

ENABLE_PARTICLE_OUTPUTS YES   # enables halo particle outputs
MIN_OUTPUT_HALO_SIZE 10000	# min size of halo to output particles from
MIN_FRACTIONAL_OUTPUT_HALO 5  # min size of halo to do fractional output
OUTPUT_PARTICLE_FRACTION 0.01 # percentage of particles to output
OUTPUT_PARTICLE_TAGS YES  	# generate a particle tags file
OUTPUT_ALL_PARTICLES NO   	# output all particles
~~~~


## Katrin’s tips for faster run:
 
* If you run this on CPUs this will take forever at step 499 (which is z=0) because the center finder is very slow. You should run it on GPUs with thrust if possible. 
 
* Another way to speed it up for now is to change the center finder option: 
~~~~
USE_MBP_FINDER  	YES	# run the MBP algorithm for FOF centers
USE_MCP_FINDER  	NO	# run the MCP algorithm for FOF centers
USE_HIST_FINDER 	NO	# run the histogram MCP algorithm for FOF centers
~~~~
could be made:
~~~~
USE_MBP_FINDER  	NO    # run the MBP algorithm for FOF centers
USE_MCP_FINDER  	NO	# run the MCP algorithm for FOF centers
USE_HIST_FINDER 	YES    # run the histogram MCP algorithm for FOF centers
This will not give you accurate centers but speed things up for testing.
~~~~
 
* Another speed up option:
MAX_FOR_CENTER_FINDING 10000 # run center finder for halos with this maximum number of particles

**Again, this will not be everything, so only for testing useful**
 
* Also, for a test, do you have a much earlier time step?