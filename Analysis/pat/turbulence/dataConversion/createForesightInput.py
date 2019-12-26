import sys
import os
import shutil 
import fileinput

#folder = "/home/pascal/projects/VizAly-Foresight/inputs/turbulence"
#startTimestep = 50
#numTimesteps = 50
#originFile = "gda-turbulence_template.json"

folder = sys.argv[1]
startTimestep = int(sys.argv[2])
numTimesteps = int(sys.argv[3])
originFile = sys.argv[4]


scalars = ["vx", "vy", "vz"]


source = folder + "/" + originFile
for i in range(numTimesteps):
	for scalar in scalars:
		destination = folder + "/" + "gda-ts" + str(startTimestep+i) + "-" + scalar + ".json"
		dest = shutil.copyfile(source, destination)


		with fileinput.FileInput(destination, inplace=True) as file:
			for line in file:
				print( line.replace("v*", scalar), end='')

		with fileinput.FileInput(destination, inplace=True) as file:
			for line in file:
				print( line.replace("t*", str(startTimestep+i)), end='')


# python createForesightInput.py /home/pascal/projects/VizAly-Foresight/inputs/turbulence 50 50 gda-turbulence_template.json