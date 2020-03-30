#! /usr/bin/env python
"""
To run executable do:

python3 -m Workflow.cfdns.workflow --input-file inputs/cfdns/ml_turbulence.json --reduction
python3 -m Workflow.cfdns.workflow --input-file inputs/cfdns/ml_turbulence_validation.json --analysis
python3 -m Workflow.cfdns.workflow --input-file inputs/cfdns/ml_turbulence_validation.json --vis
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
		
		#self.add_cbench_job()
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

			#self.add_job(analysis_job,dependencies="single", filter="postprocess")
			self.add_job(analysis_job)


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


		## Create Average metrics
		utils.average_csv(self.json_data["project-home"] + self.json_data["wflow-path"] + "/reduction/cbench/metrics", 
							self.json_data["input"]["timesteps"][1], output_filename),
							self.json_data["project-home"] + self.json_data["wflow-path"] + "/reduction/cbench/metrics_avg.csv"
  
  
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
								execute_dir="cinema",
								executable=self.json_data["visualize"]["cinema"]["path"], 
								configurations=configurations,
								environment=environment)

		# make dependent on CBench job and add to workflow
		self.add_job(cinema_job, dependencies="type", filter="plot")
		"""



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