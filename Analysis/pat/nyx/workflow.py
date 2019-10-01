#! /usr/bin/env python
"""
To run executable do:

python3 -m pat.nyx.workflow --input-file ../inputs/tests/test-nyx-foresight.json
python3 -m pat.nyx.workflow --input-file ../inputs/tests/test-nyx-analysis-cinema.json --analysis-cinema
python3 -m pat.nyx.workflow --input-file ../inputs/tests/test-nyx-cinema.json --cinema
python3 -m pat.nyx.workflow --input-file ../inputs/tests/test-nyx-foresight.json --cbench
"""

import sys
import os
import argparse

from pat.utils import workflow as workflow
from pat.utils import job as j
from pat.utils import file_utilities as futils
from pat.utils import plot_utilities as putils


class NYXWorkflow(workflow.Workflow):

	def __init__(self, name, json_data, workflow_dir=""):
		# Write the new analysis file
		super().__init__(name, json_data, workflow_dir)



	# Re-write the json data to include the analysis; ["pat"]["analysis"]
	def create_analysis_input(self):	

		if "analysis-results" in self.json_data["input"]:
			analysis_path = self.json_data["input"]["analysis-results"]
		else:
			analysis_path = self.json_data["project-home"] +  self.json_data['wflow-path']

		# Remove all entries if any
		self.json_data['pat']['analysis'].clear()


		# Add analysis entries
		for ana in self.json_data['pat']['analysis-tool']['analytics']:
			for item in ana['type']:
				json_item = {
					"title" : ana['name'] + "_" + item,
					"files" : []
				}

				for inputItem in self.json_data['pat']['input-files']:
					input_item = { 
						'name' : inputItem["output-prefix"],
						'path' : analysis_path + "/" + ana['name'] + "/" + inputItem['output-prefix'] + item + ana['postfix']
					}

					json_item['files'].append(input_item)

				self.json_data['pat']['analysis'].append(json_item)

	

	# Create the analysis job; run gimlet
	def add_analysis_jobs(self):

		# get CBench job which is parent to all jobs in this function
		has_parents = False
		if ( len(self.jobs) > 0):
			parent_job = self.jobs[0]
			has_parents = True


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
				if (has_parents):
					gimlet_job.add_parents(parent_job)

				self.add_job(gimlet_job)




	# Create plots for cinema + cinema database
	def add_cinema_plotting_jobs(self):

		self.create_analysis_input()

		if "evn_path" in self.json_data["cinema-plots"]:
			environment = self.json_data["foresight-home"] + self.json_data["cinema-plots"]["evn_path"]
		else:
			environment = None

		if "configuration" in self.json_data["cinema-plots"]:
			configurations = list( sum( self.json_data["cinema-plots"]["configuration"].items(), () ) )
		else:
			configurations = None

		arg1 = self.json_data["project-home"] +  self.json_data['wflow-path'] + "/wflow.json"
		plot_path = self.json_data['project-home'] + self.json_data['wflow-path'] + "/plots"

		cinema_job = j.Job(name="cinema_",
			execute_dir="cinema",
			executable="python -m pat.nyx.cinema", 
			arguments=[ "--input-file", arg1 ],
			configurations=configurations,
			environment=environment )
		cinema_job.add_command("mkdir " + plot_path)
		cinema_job.add_command("cd " + self.json_data["foresight-home"] + "/Analysis/")
		cinema_job.add_command("rm -rf " + self.json_data["project-home"] + "/" + self.json_data["wflow-path"] + "/" + "/cbench/" + self.json_data['cbench']['output']['output-decompressed-location'])

		# make dependent on CBench job and add to workflow
		if ( len(self.jobs) > 0):
			#print(self.jobs)
			for job in self.jobs:
				cinema_job.add_parents(job)
			self.add_job(cinema_job)
		else:
			self.add_job(cinema_job)



def main():
	# Parse Input
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("--input-file")
	parser.add_argument("--analysis-cinema", action="store_true", help="run analysis and cinema job only")
	parser.add_argument("--cbench",          action="store_true", help="run cbench only")
	parser.add_argument("--cinema",          action="store_true", help="run cinema only")
	parser.add_argument("--preview", 		 action="store_true", help="preview the job, create scripts, ... but won't run")
	opts = parser.parse_args()


	# Read input JSON file
	wflow_data = futils.read_json(opts.input_file)


	# create Workflow instance
	wflow_dir = wflow_data["project-home"] + wflow_data['wflow-path']
	wflow = NYXWorkflow("wflow", wflow_data, workflow_dir=wflow_dir)


	# add jobs to workflow
	if opts.cinema:
		print("Run cinema only")
		wflow.add_cinema_plotting_jobs()
	elif opts.cbench:
		print("Run cbench only")
		wflow.add_cbench_job()
	elif opts.analysis_cinema:
		print("Run analysis + cinema")
		wflow.add_analysis_jobs()
		wflow.add_cinema_plotting_jobs()
	else:	
		print("Run full: CBench, analysis, and Cinema")
		wflow.add_cbench_job()
		wflow.add_analysis_jobs()
		wflow.add_cinema_plotting_jobs()


	# submit workflow
	if opts.preview:
		wflow.write_submit()
	else:
		wflow.write_submit()
		wflow.submit()



if __name__== "__main__":
  	main()