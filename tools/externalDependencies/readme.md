The folder contains a scripts to build external dependencies for compressors and HDF5. Currently, there are scripts for:
 - BLOSC
 - SZ
 - BigCrunch
 - HDF5

## Adding new build scripts
 - All scripts should be prefixed with "build_" and end in ".sh"
   - if this is not the case, the script won't be executed!!!

 - Use build_BLOSC.sh as an example script
