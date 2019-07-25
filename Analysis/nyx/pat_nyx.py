#! /usr/bin/env python
"""
To run executable do:
python pat_nyx.py --input-file ../inputs/nyx/NYX_wflow.json
python pat_nyx.py --input-file ../inputs/nyx/NYX_wflow.json --submit
python pat_nyx.py --input-file ../inputs/nyx/NYX_wflow_darwin.json --submit
"""

import argparse
import os
from ..pat import file_utilities as futils
from ..pat import plot_utilities as putils
from ..pat import workflow
from ..pat import Job as j


class NYXWorkflow(workflow.Workflow):

	# Re-write the json data to include the analysis
	def add_analysis_input(self):		
		# Create analysis settings
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



	# Create the analysis job
	def add_analysis_jobs(self):

		# get CBench job which is parent to all jobs in this function
		cbench_job = self.jobs[0]

		# create job to run sim_stat and lya
		for analysis in self.json_data["pat"]["analysis-tool"]["analytics"]:

			if "evn_path" in self.json_data["pat"]:
				environment = self.json_data["foresight-home"] + self.json_data["pat"]["evn_path"]
			else:
				environment = None

			if "configuration" in analysis:
				configurations = list( sum( analysis["configuration"].items(), () ) )
			else:
				configurations = None


			for item in self.json_data["pat"]["input-files"]:
				print("Creating analysis jobs for {} on {}".format(analysis, item))

				#execute_dir=self.json_data["project-home"] + "/" + analysis["name"],
				# create job for sim_stats
				gimlet_job = j.Job(name="{}_{}".format(item["output-prefix"], analysis["name"]),
										execute_dir=analysis["name"],
										executable=analysis["path"], 
										arguments=[ item["path"], item["output-prefix"] ],
										configurations=configurations,
										environment=environment )

				if "command" in analysis:
					gimlet_job.add_command(analysis["command"])


				# make dependent on CBench job and add to workflow
				gimlet_job.add_parents(cbench_job)
				self.add_job(gimlet_job)

		self.add_analysis_input()



	# Create plots
	def add_cinema_plotting_jobs(self):
		if "evn_path" in self.json_data["cinema-plots"]:
			environment = self.json_data["foresight-home"] + self.json_data["cinema-plots"]["evn_path"]
		else:
			environment = None

		if "configuration" in self.json_data["cinema-plots"]:
			configurations = list( sum( self.json_data["cinema-plots"]["configuration"].items(), () ) )
		else:
			configurations = None

		arg1 = self.json_data["project-home"] +  self.json_data['wflow-path'] + "/cbench/wflow.json"
		plot_path = self.json_data['project-home'] + self.json_data['wflow-path'] + "/plots"

		cinema_job = j.Job(name="cinema_",
			execute_dir="cinema",
			executable="python " + self.json_data["foresight-home"] + "/Analysis/" + "pat_nyx_cinema.py", 
			arguments=[ "--input-file", arg1 ],
			configurations=configurations,
			environment=environment )
		cinema_job.add_command("mkdir " + plot_path)

		# make dependent on CBench job and add to workflow
		print(self.jobs)
		for job in self.jobs:
			cinema_job.add_parents(job)
		self.add_job(cinema_job)



# Parse Input
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--input-file")
parser.add_argument("--submit", action="store_true")
opts = parser.parse_args()

# Read input JSON file
wflow_data = futils.read_json(opts.input_file)

# create Workflow instance
wflow_dir = wflow_data["project-home"] + wflow_data['wflow-path']
wflow = NYXWorkflow("wflow", wflow_data, workflow_dir=wflow_dir)

# add jobs to workflow
wflow.add_cbench_job()
wflow.add_analysis_jobs()
wflow.add_cinema_plotting_jobs()


# write submission script
wflow.write_submit()


# submit workflow
if opts.submit:
	wflow.submit()

