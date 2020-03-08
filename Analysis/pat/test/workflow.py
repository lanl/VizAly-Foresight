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
python3 -m pat.nyx.workflow --input-file ../inputs/hacc/<>.json --cbench
"""

import argparse
import os

from pat.utils import workflow as workflow
from pat.utils import job as j
from pat.utils import utilities as utils
from pat.utils import file_utilities as futils
from pat.utils import plot_utilities as putils


class NYXWorkflow(workflow.Workflow):

	def add_data_reduction_jobs(self):
		self.add_cbench_job()
		self.fill_reduc_output_presets()


	# Re-write the json data to include the analysis
	def add_analysis_input(self):		
		""" Create analysis settings """
		analyis_path = self.json_data["project-home"] +  self.json_data['wflow-path']
		

		for ana in self.json_data['pat']['analysis-tool']['analytics']:
			for item in ana['type']:
				json_item = {
					"title" : ana['name'] + "_" + item,
					"files" : []
				}

				for inputItem in self.json_data['pat']['input-files']:
					input_item = { 
						'name' : inputItem["output-prefix"],
						'path' : analyis_path + "/" + ana['name'] + "/" + inputItem['output-prefix'] + item + ana['postfix']
					}

					json_item['files'].append(input_item)

				self.json_data['pat']['analysis'].append(json_item)


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

					if param == "%postfix%":
						current_params[i] = "_" + item["output-prefix"]


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

		#self.add_analysis_input()


	def add_cinema_plotting_jobs(self):
		""" Create the plots and draw a cinema database """

		# sources the modules to be loaded on that cluster
		if "evn_path" in self.json_data["cinema-plots"]:
			environment = self.json_data["foresight-home"] + self.json_data["cinema-plots"]["evn_path"]
		else:
			environment = None

		# pull the cluster parameters with which to launch job
		if "configuration" in self.json_data["cinema-plots"]:
			configurations = list( sum( self.json_data["cinema-plots"]["configuration"].items(), () ) )
		else:
			configurations = None


		if "commands" in self.json_data["cinema-plots"]:
			commands = self.json_data["cinema-plots"]["commands"]
		else:
			commands = None

		if "slurm_commands" in self.json_data["cinema-plots"]:
			slurm_commands = self.json_data["cinema-plots"]["slurm_commands"]
		else:
			slurm_commands = None

		# Create the job script
		arg1 = self.json_data["project-home"] +  self.json_data['wflow-path'] + "/cbench/wflow.json"
		plot_path = self.json_data['project-home'] + self.json_data['wflow-path'] + "/plots"

		cinema_job = j.Job(name="cinema_",
			execute_dir="cinema",
			executable="python " + self.json_data["foresight-home"] + "/Analysis/" + "pat_nyx_cinema.py", 
			arguments=[ "--input-file", arg1 ],
			configurations=configurations,
			environment=environment,
			commands=commands,
			slurm_commands=slurm_commands)
		
		cinema_job.add_command("mkdir " + plot_path)

		# make dependent on analysis job and add to workflow
		print(self.jobs)
		for job in self.jobs:
			cinema_job.add_parents(job)
		self.add_job(cinema_job)



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