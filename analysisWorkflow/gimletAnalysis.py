#! /usr/bin/env python
import os, sys

input_file = sys.argv[1]
output_file = sys.argv[2]
prefix = sys.argv[3]

def createCommand(hdfpath, input_file, output_file, field):
	command = hdf5_path + "h5copy" + " -v -i \"" + input_file + "\" -o \"" + output_file + "\" -s \"" + field + "\" -d \"" + field + "\""
	return command

if __name__ == "__main__":
	hdf5_path = "/projects/exasky/VizAly-CBench/ExternalDependencies/hdf5/install/bin/"

	# Copy missing attributes
	command = createCommand(hdf5_path, input_file, output_file, "/domain")
	os.system(command)

	command = createCommand(hdf5_path, input_file, output_file, "/universe")
	os.system(command)


	# setup environment
	os.system("source /projects/exasky/HACC.darwin_setup")
	os.system("cd /projects/exasky/gimlet2")

	gimlet_path = "/projects/exasky/gimlet2/apps/sim_stats/sim_stats.ex"
	cmd = gimlet_path + " " + output_file + " " + prefix
	os.system(cmd)

"""
python runAnalysis.py /projects/exasky/data/NYX/highz/512/NVB_C009_l10n512_S12345T692_z42.hdf5 /projects/exasky/VizAly-CBench/build/__SZ_1853541582__NVB_C009_l10n512_S12345T692_z42.hdf5 _test_SZ_
fields: /domain /universe
"""
