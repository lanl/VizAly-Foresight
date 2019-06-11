#!/usr/bin/python

import sys, os, json, csv

class Job:
	cluster = ""
	run_script = ""
	run_script_name = ""


	def __init__(self, cluster_name):
		self.cluter = cluster_name


	def create_run_script(self, num_nodes, num_ranks_per_nodes, 
					run_partition, load_module, run_folder, 
					exe, config, timestep, params, inputData):

		if (self.cluster == "Darwin"):
			self.run_script = self.createDarwinRunScript(num_nodes, num_ranks_per_nodes, 
														run_partition, load_module, run_folder, 
														exe, config, timestep, params, inputData)



	def create_darwin_run_script(self, num_nodes, num_ranks_per_nodes, run_partition, 
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


	def run_script(self, jobscript_name):
		self.run_script_name = jobscript_name
		run_script_file = open(jobscript_name, "w") 
		run_script_file.write(self.run_script)
		run_script_file.close()
		OS.command("sbatch " + self.run_script_name)
