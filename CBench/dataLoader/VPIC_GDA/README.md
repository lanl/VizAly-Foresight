The GDA reader is for VPIC binary data, which composed of two files: a binary file (with a .gda extension) with float data and a .info file describing the binary data with the following layout:

offset (int)
x dimensions
y dimensions
z dimensions
offset 
real x dimensions
real y dimensions
real z dimensions

An example is below:
0
511.0
511.0
511.0
0
300.0
194.0
194.0

This loader uses MPIIO to read and write binary data in parallel.