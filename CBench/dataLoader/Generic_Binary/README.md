The file format to include in the input JSON file; "filetype": "RAW"

This is a generic loader for binary data.
- the binary data is has .raw extension.
- the info file has a .info extension describing the binary data. The layout is as follows:

path to raw file
x-dimensions y-dimensions z-dimensions
scalar_name offset type
   .		  .		.
   .		  .		.

e.g.
Path to raw file
128 128 128
vx 0 float
vy 8388608 float
vx 16777216 float

Implmentation details:
 - This uses MPIIO to read and write binary data in parallel
 - multiple timesteps supported
