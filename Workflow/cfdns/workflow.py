#! /usr/bin/env python
"""
To run executable do:

python3 -m Workflow.cfdns.workflow --input-file inputs/cfdns/ml.json --reduction
"""

import sys
import os
import argparse

import Workflow.draw.workflow as workflow
import Workflow.draw.job as j
import Workflow.draw.utils as utils


class CFDNSWorkflow(workflow.Workflow):

	def __init__(self, name, json_data, workflow_dir=""):
		# Write the new analysis file
		super().__init__(name, json_data, workflow_dir)

	def preprocess_cbench(self):
		for script in self.json_data['data-reduction']['cbench-pre-process']:
      
			# Create configuration
			if "configuration" in script:
				configurations = list(sum(script["configuration"].items(), ()))
			else:
				configurations = None
    
			if "evn_path" in script["evn_path"]:
				environment =  script["evn_path"]
			else:
				environment = None

   
			# Create Job
			thisjob = j.Job(name="preprocess",
				execute_dir=script['exec-dir'],
				executable=script['path'], 
				arguments=script['params'],
				configurations=configurations,
				environment=environment)
			
			self.add_job(thisjob)


	def postprocess_cbench(self):
		for script in self.json_data['data-reduction']['cbench-post-process']:

			# Create configuration
			if "configuration" in script:
				configurations = list(sum(script["configuration"].items(), ()))
			else:
				configurations = None
    
			if "evn_path" in script["evn_path"]:
				environment =  script["evn_path"]
			else:
				environment = None

			print("script['params']:", script['params'])

			# Create Job
			thisjob = j.Job(name="postprocess",
				execute_dir=script['exec-dir'],
				executable=script['path'], 
				arguments=script['params'],
				configurations=configurations,
				environment=environment)
			
			# Add it to queue
			self.add_job(thisjob, dependencies="single", filter="cbench")
		
		
	# Run cbench job
	def add_data_reduction_jobs(self):
		#self.preprocess_cbench()
		#self.add_cbench_job(dependency="single", filters="preprocess")
		self.add_cbench_job()
		self.postprocess_cbench()


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