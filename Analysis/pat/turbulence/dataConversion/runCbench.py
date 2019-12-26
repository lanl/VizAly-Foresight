import os
import sys
import h5py
import numpy as np

foresight_path = "/home/pascal/projects/VizAly-Foresight/"
scalars = ["vx", "vy", "vz"]
startTimestep = 50
numTimesteps = 50

cmd = "pushd " + foresight_path + "build"
os.system(cmd)
os.system("pwd")

for ts in range(numTimesteps):
	for scalar in scalars:
		cmd = "export LD_LIBRARY_PATH=/home/pascal/projects/VizAly-Foresight/ExternalDependencies/SZ/install/lib:$LD_LIBRARY_PATH"
		print(cmd)
		os.system(cmd)
		cmd = "mpirun -np 8 " + foresight_path + "build/CBench " + foresight_path + "inputs/turbulence/gda-ts" + str(startTimestep+ts) + "-" + scalar + ".json"
		print(cmd)
		os.system(cmd)

os.system("popd")
