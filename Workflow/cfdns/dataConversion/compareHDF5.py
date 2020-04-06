import os
import sys
import h5py
import numpy as np


# Gather inputs
filename1 = sys.argv[1]
filename2 = sys.argv[2]

snapshots = 50
numScalars = 5
numUsedScalars = 3
dimX = 128
dimY = 128
dimZ = 128

f1 = h5py.File(filename1,'r')
fields1 = f1['fields']
#fields1 = f1['dns']
data1 = fields1[:snapshots,:,:,:,:numUsedScalars]


f2 = h5py.File(filename2,'r')
fields2 = f2['fields']
data2 = fields2[:snapshots,:,:,:,:numUsedScalars]


for ts in range(snapshots):
	print("ts: ", ts)
	for s in range(numUsedScalars):
		print("scalar: ", s)
		for z in range(dimZ):
			for y in range(dimY):
				for x in range(dimX):
					if data1[ts][x][y][z][s] != data2[ts][x][y][z][s]:
						print("Different at ", x ,", ", y ,", ", z )
						exit(0)
