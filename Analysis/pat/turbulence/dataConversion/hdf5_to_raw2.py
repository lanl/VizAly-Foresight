import os
import sys
import h5py
import numpy as np


def write_info_file(filename, dimsx, dimsy, dimsz, scalars, dataTypes):
    file = open(filename,"w")
    file.write(filename + "\n")
    file.write(str(dimsx) + "\n")
    file.write(str(dimsy) + "\n")
    file.write(str(dimsz) + "\n")
    cummulativeOffset = 0

    for i in range(len(scalars)):
        file.write(scalars[i] + " double " + str(cummulativeOffset) + "\n")
        cummulativeOffset = cummulativeOffset + dimsx * dimsy * dimsz * 8
    file.close()


def hdf5_to_raw(filename, outputFolder, scalars, timestep):
	# Reads file
	f = h5py.File(filename,'r')
	fields = f['fields']
	snapshots = 100
	data = fields[:snapshots,:,:,:,:]

	# Select timestep
	dataTimestep1 = data[timestep]

	dimx = dataTimestep1.shape[0]
	dimy = dataTimestep1.shape[1]
	dimz = dataTimestep1.shape[2]
	num_scalars = dataTimestep1.shape[3]

	print("num_scalars: ", num_scalars)

	# Select scalars
	dataTypes = []

	num_elements = dimx * dimy * dimz * num_scalars
	readInDataset = np.empty([num_elements], dtype=np.float64)

	# Process data and output
	count = 0
	for index in range(num_scalars):
		dataTypes.append( type(dataTimestep1[0][0][0][index]) )
		print("index:", index)

		for _z in range(dimz):
			for _y in range(dimy):
				for _x in range(dimx):
					readInDataset[count] = dataTimestep1[_x][_y][_z][index]
					count = count + 1
	
 	# Do output
	
 	# write file data
	outputfileName =  "turbulence_ts_" + str(timestep) + ".raw"
	out_file = open( (outputFolder+ "/" + outputfileName), 'wb')
	readInDataset.tofile(out_file)
 
	# write file info
	write_info_file(outputFolder + "/turbulence_ts_" + str(timestep) + ".info", dimx, dimy, dimz, scalars, dataTypes)

	print("wrote out ", outputfileName)



# Gather inputs
filename = sys.argv[1]

scalar_1 = sys.argv[2]
scalar_2 = sys.argv[3]
scalar_3 = sys.argv[4]
scalar_4 = sys.argv[5]
scalar_5 = sys.argv[6]

timestep1 = int(sys.argv[7])
timestep2 = int(sys.argv[8])

outputFolder = sys.argv[9]

# Read in scalars
scalars = []
scalars.append(scalar_1)
scalars.append(scalar_2)
scalars.append(scalar_3)
scalars.append(scalar_4)
scalars.append(scalar_5)


# Create Folder
try:
	os.mkdir(outputFolder)
except OSError:
	print ("Directory %s already exists" % outputFolder)
else:
	print ("Successfully created the directory %s " % outputFolder)


for ts in range(timestep1, timestep2):
	hdf5_to_raw(filename, outputFolder, scalars, ts)


# Run as:
# python dataConversion/hdf5_to_raw2.py /projects/ml_compression/data/scalarHIT_fields100.h5 one two vx vy vz 0 50  temp
