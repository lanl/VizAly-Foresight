import sys
import os
import re
import argparse
import h5py
import numpy as np



# From: https://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
def sorted_nicely(l):
	""" Sort the given iterable in the way that humans expect.""" 
 
	convert = lambda text: int(text) if text.isdigit() else text 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(l, key = alphanum_key)


def parseList(s):
	# s should be in the format ['aa','bb','cc']
	try:
		s1 = s.replace('[','')
		s2 = s1.replace(']','')
		s3 = s2.replace('\'','')

		l = list(map(str, s3.split(',')))
		return l

	except:
		raise argparse.ArgumentTypeError("Something is wrong with the list")


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
#outputFilename = sys.argv[1]

parser = argparse.ArgumentParser()

parser.add_argument('dim_x', type=int)
parser.add_argument('dim_y', type=int)
parser.add_argument('dim_z', type=int)

parser.add_argument('dataFolder')

parser.add_argument('minTimestep', type=int)
parser.add_argument('maxTimestep', type=int)

parser.add_argument('scalar1')
parser.add_argument('scalar2')
parser.add_argument('scalar3')
parser.add_argument('scalar4')
parser.add_argument('scalar5')

parser.add_argument('prefix', type=parseList)

args = parser.parse_args()


scalars = []
scalars.append(args.scalar1)
scalars.append(args.scalar2)
scalars.append(args.scalar3)
scalars.append(args.scalar4)
scalars.append(args.scalar5)


#dims_x = int(sys.argv[1])
#dims_y = int(sys.argv[2])
#dims_z = int(sys.argv[3])

#dataFolder = sys.argv[4]

#minTimestep = int(sys.argv[5])
#maxTimestep = int(sys.argv[6])

#scalars = []
#scalars.append(sys.argv[7])
#scalars.append(sys.argv[8])
#scalars.append(sys.argv[9])
#scalars.append(sys.argv[10])
#scalars.append(sys.argv[11])

#prefix = sys.argv[12]



#files = []
#for ts in range(minTimestep, maxTimestep):
#	filaname = dataFolder + "/" + prefix + str(ts) + ".raw"
#	files.append(filaname)


#prefix = ["SZ_abs_0.01_", "BLOSC_"]
#dataFolder = "/usr/projects/ml_compression/VizAly-Foresight-CBench/build/runs/test-ml-turbulence/cbench/decompressed_files"5output_filename = ""
print("prefix:", args.prefix)
print("len", len(args.prefix))

for each_prefix in args.prefix:
	print("each_prefix:", each_prefix)
	outputFilename = each_prefix + ".h5"
 
	file_list = []
	for root, dirs, files in os.walk(args.dataFolder):
		for filename in files:
			if ( filename.startswith(each_prefix) ):
				file_list.append(args.dataFolder + "/" + filename)
	
				pos = filename.find('_ts_') + 4
				output_filename = filename[:pos]
	
	files_to_process = sorted_nicely(file_list)
 
	print("Files to process: ", files_to_process)
	raw_to_hdf5(outputFilename, args.maxTimestep-args.minTimestep, args.dim_x, args.dim_y, args.dim_z, 5, files_to_process)



# Run as:
# python dataConversion/raw_to_hdf52.py testAll50.h5 128 128 128 temp 0 50 one two vx vy vz turbulence_ts_

