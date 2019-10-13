This is a generic loader for binary data.
- the binary data is has .raw extension.
- the info file has a .info extension describing the binary data. The layout is as follows:

x-dimensions y-dimensions z-dimensions
scalar_name offset type
   .		  .		.
   .		  .		.

e.g.
128 128 128
vx 0 float
vy 8388608 float
vx 16777216 float

Both the .raw and .info have the same name. Only the extension is different.

This  uses MPIIO to read and write binary data in parallel.
