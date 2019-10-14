import os
import sys
import h5py
import numpy as np


def write_info_file(filename, dimsx, dimsy, dimsz, dataType):
	file = open(filename,"w") 
 
	file.write("0\n") 
	file.write(str(dimsx) + "\n") 
	file.write(str(dimsy) + "\n") 
	file.write(str(dimsz) + "\n") 
	file.write("0\n") 
	file.write(str(dimsx) + "\n") 
	file.write(str(dimsy) + "\n") 
	file.write(str(dimsz) + "\n") 
	file.write(dataType) 	 
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

	# Select scalars
	#scalars = ["vx", "vy", "vz"]


	readInDataset = np.empty([dimx, dimy, dimz], dtype=np.float64)

	# Process data and output
	for index in range(dataTimestep1.shape[3]):
		#myList = []
		#count = 0
		for _z in range(dimz):
			for _y in range(dimy):
				for _x in range(dimx):
					readInDataset[_x][_y][_z] = dataTimestep1[_x][_y][_z][index]
					#myList.append(dataTimestep1[_x][_y][_z][index])
					#count = count + 1


		outputfileName =  "ts_" + str(timestep) + "_" + scalars[index] + ".gda"
		out_file = open( (outputFolder+ "/" + outputfileName), 'wb')
		#out_arr = np.asarray(myList)
		

		write_info_file(outputFolder + "/ts_" + str(ts) + "_" + scalars[index] + ".info", dimx, dimy, dimz, "double")
		#out_arr.tofile(out_file)
		readInDataset.tofile(out_file)
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
# python hdf5_to_raw.py /bigData/Turbulence/scalarHIT_fields100.h5 one two vx vy vz 0 50  gdaFiles 
# python3 dataConversion/hdf5_to_raw.py /bigData/Turbulence/scalarHIT_fields100.h5 one two vx vy vz 0 50  data/gdaOriginalFiles
