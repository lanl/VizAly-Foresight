#! /usr/bin/env python
""" This script generates a workflow that runs CBench and HACC data analysis executables.

python3 -m Workflow.hacc.workflow --input-file inputs/hacc/darwin_draw_wflow.json --reduction-analysis --preview
"""

import sys
import os
import argparse
import json

import Workflow.draw.workflow as workflow
import Workflow.draw.job as j
import Workflow.draw.utils as utils



class HACCWorkflow(workflow.Workflow):

	def __init__(self, name, json_data, workflow_dir=""):
		# Write the new analysis file
		super().__init__(name, json_data, workflow_dir)


	# data reduction
	def add_data_reduction_jobs(self):	
		self.add_cbench_job()


	# Analysis: halo
	def process_params(self, json_data, args):
		new_list = []
		#print("\n\nargs:", args)
		for item in args:
			#print("item:", item)

			extracted = utils.extract_string(item,'$')
			#print("extracted:", extracted)

			for e in extracted:
				#print("e:", e)

				if e != "":
					replacement = ""
					if e == "$output-prefix$":
						replacement = json_data["output-prefix"]
					elif e == "$path$":
						replacement = json_data["path"]
					elif e == "$project-home$":
						replacement = self.json_data["project-home"]
					elif e == "$wflow-path$":
						replacement = self.json_data["wflow-path"]

					item = item.replace(e, replacement)

			new_list.append(item)
		
		#print("done process_argument - new_list:", new_list)	
		return new_list

	
	
	def replaceHalo(self, tempStr):
		tokens = tempStr.split(".")
		replacement = "." + tokens[len(tokens)-1]
		newStr = tempStr.replace(replacement,"")
		return newStr
    
	def replaceData(self, strings):
		data = strings[2]
		tokens = data.split(" ")
		y =  self.replaceHalo(tokens[1])
		newStr = strings[2].replace(tokens[1],y)
		strings[2] = newStr
		return strings


	def add_analysis_jobs(self):

		# Additional processing required for halo files
		new_halo_config_path = ""
		new_parameters_path = ""
		timestep = ""
		for ana in self.json_data["analysis"]["analytics"]:
			if ana["name"] == "halo":
				halo_analysis_path = self.json_data["project-home"] + self.json_data["wflow-path"] + "/analysis/halo/"
				os.system("mkdir -p " + halo_analysis_path)


				for param in ana["params"]:
					#print("param:", param)
					if param.startswith("--config"):
						halo_config_path = param[ len("--config"): len(param) ]
						#print("halo_config_path:", halo_config_path)
						new_halo_config_path = self.json_data["project-home"] + self.json_data["wflow-path"] + "/analysis/halo/cosmotools_config.dat"

						os.system("cp " + halo_config_path + " " + new_halo_config_path)
					
					if param.endswith("indat.params"):
						parameters_path = param
						new_parameters_path = self.json_data["project-home"] + self.json_data["wflow-path"] + "/analysis/halo/indat.params"
						os.system("cp " + parameters_path + " " + new_parameters_path)

					if param.startswith("--timesteps"):
						timestep = utils.read_file( param[ len("--timesteps")+1: len(param) ] )
		
			

		count = 0
		for output in self.json_data["data-reduction"]["output-files"]:

			# create job to run sim_stat and lya
			for analysis in self.json_data["analysis"]["analytics"]:

				configurations = utils.get_configuration( analysis )
				environment = utils.get_environment( analysis, self.json_data["foresight-home"] )
				args = self.process_params( output, analysis['params'])

				# Account for some quikeness of HACC
				if analysis["name"] == "halo":
					args = self.replaceData(args)

				
				if analysis["name"] == "halo":
					index = 0
					for param in args:
						if param.startswith("--config"):
							current_halo_config_path = self.json_data["project-home"] + self.json_data["wflow-path"] + "/analysis/halo/" + output["output-prefix"] + ".dat"
							#print("new_param_path", current_halo_config_path)
							os.system("cp " + new_halo_config_path + " " + current_halo_config_path)

							utils.replace_line_starting_with_in_file(current_halo_config_path, "BASE_OUTPUT_FILE_NAME", "BASE_OUTPUT_FILE_NAME " + output["output-prefix"] + "\n")
							utils.replace_line_starting_with_in_file(current_halo_config_path, "EXPLICIT_TIMESTEPS", "EXPLICIT_TIMESTEPS " + timestep +"\n")
							utils.replace_line_starting_with_in_file(current_halo_config_path, "ACCUMULATE_CORE_NAME", "ACCUMULATE_CORE_NAME " +  output["output-prefix"] +"\n")

							args[index] = "--config " + output["output-prefix"] + ".dat"



						if param.endswith("indat.params"):
							current_comsotools_path = self.json_data["project-home"] + self.json_data["wflow-path"] + "/analysis/halo/" + output["output-prefix"] + ".params"
							#print("current_comsotools_path", current_comsotools_path)
							os.system("cp " + new_parameters_path + " " + current_comsotools_path)

							utils.replace_line_starting_with_in_file(current_comsotools_path, "COSMOTOOLS_CONFIG", "COSMOTOOLS_CONFIG " + output["output-prefix"] + ".dat\n")

							args[index] = output["output-prefix"] + ".params"


						index = index + 1



				analysis_job = j.Job( name=analysis["name"] + str(count),
										job_type = "analysis",
										execute_dir="analysis/" + analysis["name"],
										executable=analysis["path"], 
										arguments=args,
										configurations=configurations,
										environment=environment )

				job_Dependecies = utils.get_jobDependency( analysis )
				self.add_job( analysis_job, dependencies=job_Dependecies[0], filter=job_Dependecies[1] )
			
			count = count + 1

		self.create_analysis_output()



	def create_analysis_output(self):
		self.json_data['analysis']['output-files'].clear()
	
		base_path = self.json_data['project-home'] + self.json_data['wflow-path']


		for analysis in self.json_data["analysis"]["analytics"]:

			if analysis["name"] == "halo":
				json_item = analysis["name"] : []
				self.json_data['analysis']['output-files'].append(json_item)


				# for output in self.json_data["data-reduction"]["output-files"]:
				# 	# Create output					
				# 	json_item = {
				# 		"output-prefix" : output["output-prefix"],
				# 		"path" : base_path + "/analysis/halo/" + output["output-prefix"] + "-499.haloproperties"
				# 	}

				# 	print(json_item)
				# 	self.json_data['analysis']['output-files']['halo'].append(json_item)


			elif analysis["name"] == "spectrum":
				# json_item = {
				# 	analysis["name"] : []
				# }
				# self.json_data['analysis']['output-files'].append(json_item)


			# 	for output in self.json_data["data-reduction"]["output-files"]:
			# 		# Create output					
			# 		json_item = {
			# 			"output-prefix" : output["output-prefix"],
			# 			"path" : base_path + "/analysis/spectrum/" + output["output-prefix"] + "_.pk"
			# 		}

			# 		self.json_data['analysis']['output-files']['spectrum'].append(json_item)

			else:
				# json_item = {
				# 	analysis["name"] : []
				# }
				# self.json_data['analysis']['output-files'].append(json_item)
				
			# 	self.json_data['analysis']['output-files'].append(analysis["name"])
			# 	analysis_name = analysis["name"]

			# 	for output in self.json_data["data-reduction"]["output-files"]:
			# 		# Create output					
			# 		json_item = {
			# 			"output-prefix" : output["output-prefix"],
			# 			"path" : base_path + "/analysis/" + analysis_name + "/" + output["output-prefix"] + ".txt"
			# 		}

			# 		self.json_data['analysis']['output-files'][analysis_name].append(json_item)



	def add_vis_jobs(self):
		
		## Create Plots
		for plot in self.json_data["visualize"]["plots"]:

			# sources the modules to be loaded on that cluster
			if "evn_path" in plot:
				environment = self.json_data["foresight-home"] + plot["evn_path"]
			else:
				environment = None


			# pull the cluster parameters with which to launch job
			if "configuration" in plot:
				configurations = list( sum( plot["configuration"].items(), () ) )
			else:
				configurations = None

			plot_job = j.Job(name=plot["name"],
									job_type = "plot",
									execute_dir="plot/" + plot["name"],
									executable=plot["path"], 
									arguments=plot['params'],
									configurations=configurations,
									environment=environment )

			#self.add_job(plot_job,dependencies="type", filter="analysis")
			self.add_job(plot_job)




def main():
	# Parse Input
	opts = workflow.parse_args()
 
	# Read input JSON file
	wflow_data = utils.read_json(opts.input_file)
	if wflow_data == None:
		sys.exit(0)


	# create Workflow instance
	wflow_dir = wflow_data["project-home"] + wflow_data["wflow-path"]
	wflow = HACCWorkflow("wflow", wflow_data, workflow_dir=wflow_dir)
 
	# process intput
	wflow = workflow.process_input(wflow,opts)

	# submit workflow
	if opts.preview:
		wflow.write_submit()
	else:
		wflow.write_submit()
		wflow.submit()


if __name__== "__main__":
	main()
