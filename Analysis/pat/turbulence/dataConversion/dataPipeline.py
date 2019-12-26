import os
import sys
import h5py
import numpy as np


inputHDFData = sys.argv[1]
outputRawFolder = sys.argv[2]
templateJSONFile = sys.argv[3]
foresightInputFile = sys.argv[4]
outputHDFData = sys.argv[5]


# Step 1: Convert the data from HDF5 to binary
cmd = "python3 hdf5_to_raw.py " + inputData + " one two vx vy vz 50 100 " +  outputRawFolder
print(cmd)
#os.system(cmd)

# Step 2: Create Foresight Sweep input
cmd = "python3 createForesightInput.py " + templateJSONFile + "50 50 " + foresightInputFile
print(cmd)
#os.system(cmd)

# Step 3: Run Foresight
cmd = "python3 runCbench.py"
print(cmd)
#os.system(cmd)

# Step 4: Convert back to hdf
cmd = "python3 raw_to_hdf5.py " + outputHDFData + " 128 128 128 " + decompressed_files + "50 100 one two vx vy vz" 
print(cmd)
#os.system(cmd)