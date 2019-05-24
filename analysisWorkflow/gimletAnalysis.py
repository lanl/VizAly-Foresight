#! /usr/bin/env python

import sys, json, os


def createCommand(hdfpath, input_file, output_file, field):
	command = hdf5_path + "h5copy" + " -v -i \"" + input_file + "\" -o \"" + output_file + "\" -s \"" + field + "\" -d \"" + field + "\""
	return command



if __name__ == "__main__":
	if len(sys.argv) < 2:
		print ("Json file and title needed; e.g. python gimletAnalysis.py gimletInput.json")
		exit()

	# Open json file
	with open(sys.argv[1], "r") as read_file:
		json_data = json.load(read_file)

	
	# Find original file
	index_orig = 0
	for file in json_data["input-files"]:
		if (file["name"] != "orig"):
			break

		index_orig = index_orig + 1


	# Set some environment variables
	hdf5_path = json_data["hdf5_path"]
	os.system("source " + json_data["evn_path"])
	os.system("cd " + json_data["gimlet-home"])


	# Copy uncompressed attributes from original to compressed; needed for gimlet
	for file in json_data["input-files"]:
		if (file["name"] != "orig"):
			# Copy missing attributes
			command = createCommand(hdf5_path, json_data["input-files"][index_orig]["path"], file["path"], "/domain")
			print command
			os.system(command)

			command = createCommand(hdf5_path, json_data["input-files"][index_orig]["path"], file["path"], "/universe")
			print command
			os.system(command)

		# Run gimlet analysis
		cmd = json_data["gimlet-path"] + " " + file["path"] + " " + file["output-prefix"]
			print cmd
			os.system(cmd)


"""
python runAnalysis.py gimletInput.json
"""