#!/usr/bin/python

import sys, json, os



def splitString(filename, char):
	k = filename.rfind(char)
	return filename[:k+1], filename[k+1:]



def createDarwinRunScript(num_nodes, num_ranks_per_nodes, run_partition, 
							load_module, run_folder, exe, 
							config, timestep, params, inputData):

	run_script = "#!/bin/bash\n\n"

	run_script += "#SBATCH -N " + str(num_nodes) + "\n"
	run_script += "#SBATCH --ntasks-per-node " + str(num_ranks_per_nodes) + "\n"
	run_script += "#SBATCH -p " + run_partition + "\n\n"

	# load modules
	run_script += "source " + load_module + "\n\n"

	# go to folder 
	run_script += "cd " + run_folder + "\n\n"

	# Runc command
	run_script += "mpirun " + exe
	run_script += " --config " + config
	run_script += " --timesteps " + timestep
	run_script += " --prefix " + inputData + " " + params
	run_script += " "

	return run_script



def output_script(jobscript, run_script):
	run_script_file = open(jobscript, "w") 
	run_script_file.write(run_script)
	run_script_file.close()



if __name__ == "__main__":
	# check if the json file is here
	if len(sys.argv) < 2:
		print ("Json file needed")
		exit()

	# parse Json file
	with open(sys.argv[1], "r") as read_file:
		json_data = json.load(read_file)


	# Create Path for Halo output if it doesn't exist
	halo_output_directory = json_data["halo-output-path"]
	if not os.path.exists(halo_output_directory):
    	os.makedirs(halo_output_directory)


	# Run analysis for each input provided 
	for dataset_file in json_data["input-files"]:

		_config_filename = json_data["halo-finder-setup"]["config"]
		_params_filename = json_data["halo-finder-setup"]["params"]

		config_path, config_filename = splitString(json_data["halo-finder-setup"]["config"], '/')
		params_path, params_filename = splitString(json_data["halo-finder-setup"]["params"], '/')

		current_path = os.getcwd() + '/'

		config_filename = "__" + dataset_file["alias"] + "__" + config_filename
		params_filename = "__" + dataset_file["alias"] + "__" + params_filename


		# Create new indat.param file based on indat.param file provided
		os.system('sed \'s/^COSMOTOOLS_CONFIG.*/COSMOTOOLS_CONFIG ' + config_filename + '/\' ' + _params_filename + ' > ' + params_filename)

		# Create new cosmotools-config- file based on cosmotools-config- file provided
		os.system('sed \'s/^BASE_OUTPUT_FILE_NAME.*/BASE_OUTPUT_FILE_NAME ..\/run\/output\/analysis\/Halos\/b0168\/__' + dataset_file["alias"] + '__m001/\' ' + _config_filename + ' > temp_file')
		os.system('sed \'s/^ACCUMULATE_CORE_NAME.*/ACCUMULATE_CORE_NAME ..\/run\/output\/analysis\/Halos\/b0168\/__' + dataset_file["alias"] + '__m001/\' temp_file > ' + config_filename)

		# create run script
		darwin_run_script = createDarwinRunScript(
							json_data["job-running-params"]["num-nodes"],
							json_data["job-running-params"]["num-ranks-per-nodes"],
							json_data["job-running-params"]["partition"],
							json_data["job-running-params"]["path-to-environment-setup"],
							json_data["job-running-params"]["running-path"],
							json_data["job-running-params"]["path-hacc-exe"],
							current_path + config_filename,
							json_data["halo-finder-setup"]["timesteps"],
							current_path + params_filename,
							dataset_file["name"])

		# output script
		jobscript = "job" + "_" + dataset_file["alias"] 
		output_script(jobscript, darwin_run_script)
	
		# run script file
		os.system('sbatch ' + jobscript)

# Run as:
# python workflow_a.py workflow_input.json
