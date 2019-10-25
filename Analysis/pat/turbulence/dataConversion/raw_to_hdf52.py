import sys
import h5py
import numpy as np


def raw_to_hdf5(outputFilename, numTimesteps, dimsx, dimsy, dimsz, numScalars, inputDataFiles):
	
	# Create an empty dataset
	toWtriteDataset = np.empty([numTimesteps, dimsx, dimsy, dimsz, numScalars], dtype=np.float64)
	print("num scalars: ", numScalars)

	timeIndex = 0
	for fname in inputDataFiles:
		print("Reading in:", fname, ", ts:", timeIndex)

		readData = np.fromfile(fname, dtype=np.float64)
		num_elems = dimsx * dimsy * dimsz * numScalars
		resahpedData = readData.reshape(num_elems)
		
		index = 0
		for s in range(numScalars):
			for z in range(dimsz):
				for y in range(dimsy):
					for x in range(dimsx):
						toWtriteDataset[timeIndex][x][y][z][s] = resahpedData[index]
						index = index + 1

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

scalars = []
scalars.append(sys.argv[8])
scalars.append(sys.argv[9])
scalars.append(sys.argv[10])
scalars.append(sys.argv[11])
scalars.append(sys.argv[12])

prefix = sys.argv[13]


files = []
for ts in range(minTimestep, maxTimestep):
	filaname = dataFolder + "/" + "turbulence_ts_" + str(ts) + ".raw"
	files.append(filaname)

raw_to_hdf5(outputFilename, maxTimestep-minTimestep, dims_x, dims_y, dims_z, 5, files)



# Run as:
# python dataConversion/raw_to_hdf52.py testAll50.h5 128 128 128 temp 0 50 one two vx vy vz turbulence_ts_

