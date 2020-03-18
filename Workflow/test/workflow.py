#! /usr/bin/env python
"""
pat_nyx fun the foresight pipeline for nyx analysis. In other words, it does:
- compression
- data analysis
- plot graphs
- create a cinema data

To run:
python3 -m Workflow.test.workflow --input-file inputs/test/test_nyx_darwin.json --vis
"""

import argparse
import os
import sys

import Workflow.draw.workflow as workflow
import Workflow.draw.job as j
import Workflow.draw.utils as utils


class NYXWorkflow(workflow.Workflow):

	def add_data_reduction_jobs(self):
		self.add_cbench_job()
		#self.fill_reduc_output_presets()


	# Re-write the json data to include the analysis; ["pat"]["analysis"]
	def fill_analysis_results(self):	

		analysis_path = self.json_data["project-home"] +  self.json_data['wflow-path'] + "analysis/"

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
						current_params[i] = item["output-prefix"]


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

		"""
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

			# Parameters to run the job with
			params = utils.get_list_from_json_array(plot, "params")


			# create job for sim_stats
			plot_job = j.Job(name=plot["name"],
								job_type = "plot",
								execute_dir=plot["name"] +"/",
								executable=plot["path"], 
								arguments=params,
								configurations=configurations,
								environment=environment)
 
			# make dependent on CBench job and add to workflow
			#self.add_job(plot_job, dependencies="type", filter="analysis")
			self.add_job(plot_job)
		"""
  
		"""
  "cinema" : 
		{
			"name" : "nyx_cinema",
			"path" : "python -m $foresight-home$/Workflow/test/cinema.py",
			"params" : ["/projects/exasky/pascal-projects/VizAly-Foresight/inputs/test/test_nyx_darwin.json"],
			"evn_path": "scripts/VizAly-CBench.bash.darwin",
			"configuration": 
			{
				"partition": "general",
				"nodes": 1,
				"ntasks-per-node": 1
			}
		}
		"""
  
		## Create Cinema DB
		if "evn_path" in self.json_data["visualize"]["cinema"]:
			environment = self.json_data["visualize"]["cinema"]["evn_path"] 
		else:
			environment = None

		# pull the cluster parameters with which to launch job
		if "configuration" in self.json_data["visualize"]["cinema"]:
			configurations = list( sum( self.json_data["visualize"]["cinema"]["configuration"].items(), () ) )
		else:
			configurations = None
  
  
		# create job for sim_stats
		cinema_job = j.Job(name=self.json_data["visualize"]["cinema"]["name"],
								job_type = "cinema",
								execute_dir="",
								executable=self.json_data["visualize"]["cinema"]["path"], 
								configurations=configurations,
								environment=environment)

		# make dependent on CBench job and add to workflow
		#self.add_job(cinema_job, dependencies="type", filter="plot")
		self.add_job(cinema_job)



def main():
	# Parse Input
	opts = workflow.parse_args()
 
	# Read input JSON file
	wflow_data = utils.read_json(opts.input_file)
	if wflow_data == None:
		sys.exit(0)


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