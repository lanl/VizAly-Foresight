import sys
import h5py
import numpy as np


def raw_to_hdf5(outputFilename, numTimesteps, dimsx, dimsy, dimsz, numScalars, inputDataFiles):
	
	# Create an empty dataset
	toWtriteDataset = np.empty([numTimesteps, dimsx, dimsy, dimsz, numScalars], dtype=np.float64)

	timeIndex = 0
	fieldIndex = 0
	for fname in inputDataFiles:
		print("Reading in ", fname, " ts ", timeIndex, " field index: ", fieldIndex)

		readData = np.fromfile(fname, dtype=np.float64)
		resahpedData = readData.reshape(dimsx, dimsy, dimsz)
		
		for z in range(dimsz):
			for y in range(dimsy):
				for x in range(dimsx):
					toWtriteDataset[timeIndex][x][y][z][fieldIndex] = resahpedData[x][y][z]

		fieldIndex = fieldIndex + 1
		if fieldIndex >= numScalars:
			fieldIndex = 0
			timeIndex = timeIndex + 1

	f = h5py.File(outputFilename,'w')
	f.create_dataset("fields", data=toWtriteDataset)


# Gather inputs
outputFilename = sys.argv[1]

dims_x = int(sys.argv[2])
dims_y = int(sys.argv[3])
dims_z = int(sys.argv[4])

dataFolder = sys.argv[5]

minTimestep = int(sys.argv[6])
maxTimestep = int(sys.argv[7])

scalar_one = sys.argv[8]
scalar_two = sys.argv[9]
scalar_x = sys.argv[10]
scalar_y = sys.argv[11]
scalar_z = sys.argv[12]

prefix = sys.argv[13]

files = []
for ts in range(minTimestep, maxTimestep):
	filaname = dataFolder + "/ts_" + str(ts) + "_" + scalar_one + ".gda"
	files.append(filaname)

	filaname = dataFolder + "/ts_" + str(ts) + "_" + scalar_two + ".gda"
	files.append(filaname)

	filaname = dataFolder + "/" + prefix + str(ts) + "_" + scalar_x + ".gda"
	files.append(filaname)

	filaname = dataFolder + "/" + prefix + str(ts) + "_" + scalar_y + ".gda"
	files.append(filaname)

	filaname = dataFolder + "/" + prefix + str(ts) + "_" + scalar_z + ".gda"
	files.append(filaname)


raw_to_hdf5(outputFilename, maxTimestep-minTimestep, dims_x, dims_y, dims_z, 5, files)



# Run as:
# python raw_to_hdf5.py testAll50.h5 128 128 128 decompressed_files 0 50 one two vx vy vz ts_

