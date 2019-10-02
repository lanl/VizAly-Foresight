NYX uses the HDF5 file format. It contains a number of grids. 

This loader will allow you to read and write any binary HDF in parallel.

**Note:** Parallel HDF5 needs a parallel file system for writing out data. This code will **NOT** work in Parallel on Darwin which does not have a parallel file system!