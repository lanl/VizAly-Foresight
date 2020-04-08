#! /usr/bin/env python
""" This script generates a workflow that runs CBench and HACC data analysis executables.

python3 -m pat.hacc.workflow --input-file ../inputs/hacc/<>.json-foresight.json
python3 -m pat.hacc.workflow --input-file ../inputs/hacc/<>.json --analysis-cinema
python3 -m pat.hacc.workflow --input-file ../inputs/hacc/<>.json --cinema
python3 -m pat.hacc.workflow --input-file ../inputs/hacc/<>.json --cbench
"""

import sys
import os
import argparse

from pat.utils import workflow as workflow
from pat.utils import job as j
from pat.utils import file_utilities as futils


class HACCWorkflow(workflow.Workflow):

	def __init__(self, name, json_data, workflow_dir=""):
		# Write the new analysis file
		super().__init__(name, json_data, workflow_dir)


	# data reduction
	def add_data_reduction_jobs(self):	
		self.add_cbench_job()


	# Analysis: halo, 
	def add_analysis_jobs(self):
		# create job to run sim_stat and lya
		for analysis in self.json_data["analysis"]["analytics"]:

			if "evn_path" in analysis:
				environment = self.json_data["foresight-home"] + analysis["evn_path"]
			else:
				environment = None

			if "configuration" in analysis:
				configurations = list( sum( analysis["configuration"].items(), () ) )
			else:
				configurations = None


			analysis_job = j.Job(name=analysis["name"],
									job_type = "analysis",
									execute_dir="analysis/" + analysis["name"],
									executable=analysis["path"], 
									arguments=analysis['params'],
									configurations=configurations,
									environment=environment )

			self.add_job(analysis_job,dependencies="single", filter="cbench")
			#self.add_job(analysis_job)


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
	wflow = CFDNSWorkflow("wflow", wflow_data, workflow_dir=wflow_dir)
 
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
