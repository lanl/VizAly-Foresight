#! /usr/bin/env python
"""
pat_nyx fun the foresight pipeline for nyx analysis. In other words, it does:
- run cbench
- run the analysis tool
- create a cinema databse of the result

To run:
python pat_nyx.py --input-file <absolute path to input file> <--submit>

--submit: causes the pipeline to be launched. Without it, you can see the scripts generated for debugging.

e.g.
python -m tests.workflow --input-file
"""

import argparse
import os

from draw.workflow import workflow as workflow
from draw.job import job as j
from draw.utilities import utilities as utils
from draw.plot_utilities import plot_utilities as putils


class NYXWorkflow(workflow.Workflow):

	def add_data_reduction_jobs(self):
		self.add_cbench_job()
		#self.fill_reduc_output_presets()


	# Re-write the json data to include the analysis; ["pat"]["analysis"]
	def fill_analysis_results(self):	

		analysis_path = self.json_data["project-home"] +  self.json_data['wflow-path']

		# Remove all entries if any
		self.json_data['analysis']['output-files'].clear()

		# Add analysis entries
		for ana in self.json_data['analysis']['output-files-filter']:
			for item in ana['fields']:
				json_item = {
					"title" : ana['name'] + "_" + item,
					"files" : []
				}

				for inputItem in self.json_data['data-reduction']['output-files']:
					input_item = { 
						'name' : inputItem["output-prefix"],
						'path' : analysis_path + "/" + ana['name'] + "/" + inputItem['output-prefix'] + item + ana['postfix']
					}

					json_item['files'].append(input_item)

				self.json_data['analysis']['output-files'].append(json_item)


	def add_analysis_jobs(self):
		""" specify the gimlet jobs to run """

		for analysis in self.json_data["analysis"]["analytics"]:

			# sources the modules to be loaded on that cluster
			if "evn_path" in analysis:
				environment = self.json_data["foresight-home"] + analysis["evn_path"]
			else:
				environment = None


			# pull the cluster parameters with which to launch job
			if "configuration" in analysis:
				configurations = list( sum( analysis["configuration"].items(), () ) )
			else:
				configurations = None


			if "commands" in analysis:
				cmds = analysis["commands"]
			else:
				cmds = []


			# Parameters to run the job with
			params = utils.get_list_from_json_array(analysis, "params")


			# Create the job script
			for item in self.json_data["data-reduction"]["output-files"]:
				print("Creating analysis jobs for {} on {}".format(analysis, item))

				# Replace %data% param by item
				current_params = params.copy()
				for i, param in enumerate(current_params):
					if param == "%data%":
						current_params[i] = item["path"]

					if param == "%prefix%":
						current_params[i] = item["output-prefix"] + "_"


				# create job for sim_stats
				analysis_job = j.Job(name="{}_{}".format(item["output-prefix"], analysis["name"]),
										job_type = "analysis",
										execute_dir="analysis/" + analysis["name"],
										executable=analysis["path"], 
										arguments=current_params,
										configurations=configurations,
										environment=environment)

				if cmds != []:
					for cmd in cmds:
						analysis_job.add_command(cmd)
 
				# make dependent on CBench job and add to workflow
				self.add_job(analysis_job,dependencies="single", filter="cbench")
		
		self.fill_analysis_results()



	def add_vis_jobs(self):
		""" Create the plots and draw a cinema database """

		## Create Plots
		for plot in self.json_data["plots"]["visualize"]:

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


			# Parameters to run the job with
			params = utils.get_list_from_json_array(analysis, "params")


			for item in self.json_data["analysis"]["output-files"]:
				print("Creating plots jobs for {} on {}".format(analysis, item))

				if plot["name"] in item["title"]:
				
					current_params = params.copy()
					for i, param in enumerate(current_params):
						if param == "%title%":
							current_params[i] = item["title"]

						if param == "%files%":
							file_list = []
							for f in item["files"]
								file_list.append([f["name"],f["path"]])
							current_params[i] = file_list

				# create job for sim_stats
				plot_job = j.Job(name=plot["name"]),
										job_type = "plot",
										execute_dir="plots/" + plot["name"],
										executable=plot["path"], 
										arguments=current_params,
										configurations=configurations,
										environment=environment)
 
				# make dependent on CBench job and add to workflow
				self.add_job(analysis_job, dependencies="type", filter="analysis")


		## Create Cinema DB





def main():
	# Parse Input
	opts = workflow.parse_args()
 
	# Read input JSON file
	wflow_data = utils.read_json(opts.input_file)

	# create Workflow instance
	wflow_dir = wflow_data["project-home"] + wflow_data["wflow-path"]
	wflow = NYXWorkflow("wflow", wflow_data, workflow_dir=wflow_dir)
 
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