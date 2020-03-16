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


	# Re-write the json data to include the analysis; ["pat"]["analysis"]
	def create_analysis_input(self):	

		if "analysis-results" in self.json_data["input"]:
			analysis_path = self.json_data["input"]["analysis-results"]
		else:
			analysis_path = self.json_data["project-home"] +  self.json_data['wflow-path']

		# Remove all entries if any
		#self.json_data['pat']['analysis'].clear()


		# Add analysis entries
		#for ana in self.json_data['pat']['analysis-tool']['analytics']:
		#	for item in ana['type']:
		#		json_item = {
		#			"title" : ana['name'] + "_" + item,
		#			"files" : []
		#		}

		#		for inputItem in self.json_data['pat']['input-files']:
		#			input_item = { 
		#				'name' : inputItem["output-prefix"],
		#				'path' : analysis_path + "/" + ana['name'] + "/" + inputItem['output-prefix'] + item + ana['postfix']
		#			}
		#
		#			json_item['files'].append(input_item)
		#
		#		self.json_data['pat']['analysis'].append(json_item)


	# Create the analysis job
	def add_analysis_jobs(self):

		# get CBench job which is parent to all jobs in this function
		#cbench_job = self.jobs[0]
		# get CBench job which is parent to all jobs in this function
		has_parents = False
		if ( len(self.jobs) > 0):
			parent_job = self.jobs[0]
			has_parents = True


		# get environment script
		environment = self.environment_from_json_data()


		# get halo finder information from configuration file
		halo_section = "halo"
		halo_config = self.configuration_from_json_data(halo_section)
		halo_exe = self.json_data["pat"]["analysis-tool"]["analytics"][halo_section]["path"]
		timesteps_path = self.json_data["pat"]["analysis-tool"]["analytics"][halo_section]["timesteps-file"]
		config_path = self.json_data["pat"]["analysis-tool"]["analytics"][halo_section]["config-file"]
		parameters_path = self.json_data["pat"]["analysis-tool"]["analytics"][halo_section]["parameters-file"]


		# get halo distribution information from configuraiton file
		halo_query_section = "halo_query"
		halo_query_config = self.configuration_from_json_data(halo_query_section)
		halo_query_exe = self.json_data["pat"]["analysis-tool"]["analytics"][halo_query_section]["path"]


		# get spectrum information from configuration file
		spectrum_section = "spectrum"
		spectrum_config = self.configuration_from_json_data(spectrum_section)
		spectrum_exe = self.json_data["pat"]["analysis-tool"]["analytics"][spectrum_section]["path"]
		spectrum_config_path = self.json_data["pat"]["analysis-tool"]["analytics"][spectrum_section]["config-file"]


		# create jobs to run analysis jobs on each output file from CBench
		for path in self.json_data["pat"]["input-files"]:
			print("Creating analysis jobs for", path)
			prefix = path["output-prefix"]
			cbench_path = path["path"]
			timestep = cbench_path.split(".")[-1]


			# write halo finder parameters file
			# specify location of parsed configuration file inside
			new_parameters_path = halo_section + "/" + prefix + "_halo_params.txt"
			new_config_path = halo_section + "/" + prefix + "_halo_config.txt"
			os.system("sed \"s/^COSMOTOOLS_CONFIG.*/COSMOTOOLS_CONFIG .\/{}/\" {} > {}".format(
									os.path.basename(new_config_path),
									parameters_path,
									self.workflow_dir + "/" + new_parameters_path))

			# write halo finder timestep file
			# specify from parsing the input file
			new_timesteps_path = halo_section + "/" + prefix + "_halo_timesteps.txt"
			tfile = open(self.workflow_dir + "/" + new_timesteps_path, 'w+')
			tfile.write(timestep)
			tfile.close()

			# write halo finder configuration file
			# specify output prefix inside
			tmp_path = "tmp.out"
			tmp_path2 = "tmp2.out"
			os.system("sed \"s/^BASE_OUTPUT_FILE_NAME.*/BASE_OUTPUT_FILE_NAME .\/{}/\" {} > {}".format(
									prefix, config_path, tmp_path))
			os.system("sed \"s/^EXPLICIT_TIMESTEPS.*/EXPLICIT_TIMESTEPS {}/\" {} > {}".format(
                                    timestep, tmp_path, tmp_path2))
			os.system("sed \"s/^ACCUMULATE_CORE_NAME.*/ACCUMULATE_CORE_NAME .\/{}/\" {} > {}".format(
									prefix, tmp_path2, self.workflow_dir + "/" + new_config_path))
			os.remove(tmp_path)


			# create job for halo finder
			halo_job = j.Job(name="{}_{}".format(prefix, halo_section),
							 execute_dir=halo_section,
							 executable=halo_exe,
							 arguments=["--config", os.path.basename(new_config_path),
										"--timesteps", os.path.basename(new_timesteps_path),
										"--prefix", cbench_path[:-len(str(timestep)) - 1],
										os.path.basename(new_parameters_path)],
							 configurations=halo_config,
							 environment=environment)


			# make dependent on CBench job and add to workflow
			#halo_job.add_parents(cbench_job)
			if (has_parents):
				halo_job.add_parents(parent_job)
			self.add_job(halo_job)


			# predict halo finder file
			halo_file = self.workflow_dir + "/" + halo_section + "/" + prefix + "-" + timestep + ".fofproperties"


			# loop over different halo SQL queries
			# create job for getting halo distribution
			sec = self.json_data["pat"]["analysis-tool"]["analytics"][halo_query_section]
			for i, _ in enumerate(sec["query"]):
				halo_query_file = self.workflow_dir + "/" + halo_query_section + "/halo_query_{}_{}.csv".format(prefix, i)
				args = ["--input-file", halo_file,
						"--output-file", halo_query_file,
						"--query", "\"{}\"".format(sec["query"][i]),
						"--xlabel", "\"{}\"".format(sec["xlabel"][i]),
						"--ylabel", "\"{}\"".format(sec["ylabel"][i])]
				if "xlim" in sec.keys():
					args += ["--xlim", " ".join(map(str, sec["xlim"][i]))]
				if "log-bins" in sec.keys():
					if sec["log-bins"][i]:
						args += ["--log-bins"]
				halo_query_job = j.Job(name="{}_{}_{}".format(prefix, halo_query_section, i),
									   execute_dir=halo_query_section,
									   executable=halo_query_exe,
									   arguments=args,
									   configurations=halo_query_config,
									   environment=environment)
				self.json_data["pat"]["analysis"].append({"output-column" : "FILE_Halo_Distribution_{}".format(i),
														  "output-prefix" : prefix, "path" : halo_query_file})   
 

				# make dependent on halo finder job
				halo_query_job.add_parents(halo_job)
				self.add_job(halo_query_job)

			# TODO: Automate 'TOPOLOGY' in /projects/exasky/HACC/run/inputs/indat.params
			# create job for power spectrum
			spectrum_job = j.Job(name="{}_{}".format(prefix, spectrum_section),
								 execute_dir=spectrum_section,
								 executable=spectrum_exe,
								 arguments=[spectrum_config_path, "-n", cbench_path, prefix + "_spectrum", timestep],
								 configurations=spectrum_config,
								 environment=environment)


			# make dependent on CBench job and add to workflow
			#spectrum_job.add_parents(cbench_job)
			if (has_parents):
				spectrum_job.add_parents(parent_job)
			self.add_job(spectrum_job)


			# predict and track output files
			for ext in ["", ".rsd.0", ".rsd.1", ".rsd.2"]:
				spectrum_file = self.workflow_dir + "/" + spectrum_section + "/{}_spectrum.pk{}".format(prefix, ext)
				self.json_data["pat"]["analysis"].append({"output-column" : "FILE_Spectrum{}".format(ext),
														  "output-prefix" : prefix, "path" : spectrum_file})


	# Create plots for cinema + cinema database
	def add_cinema_plotting_jobs(self):

		#self.create_analysis_input()

		# get environment for Cinema job
		if "evn_path" in self.json_data["cinema-plots"]:
			environment = self.json_data["foresight-home"] + "/"+ self.json_data["cinema-plots"]["evn_path"]
		else:
			environment = None

		# get configuration for Cinema job
		if "configuration" in self.json_data["cinema-plots"]:
			configurations = list(sum(self.json_data["cinema-plots"]["configuration"].items(), ()))
		else:
			configurations = None

		arg1 = self.json_data["project-home"] +  self.json_data['wflow-path'] + "/wflow.json"
		arg2 = self.json_data["project-home"] +  self.json_data['wflow-path'] + "/cinema/results.cdb"
		plot_path = self.json_data['project-home'] + self.json_data['wflow-path'] + "/plots"

		# get executable for Cinema job
		cinema_exe = self.json_data["cinema-plots"]["path"]

		# create job for Cinema
		cinema_job = j.Job(name="cinema",
						   execute_dir="cinema",
						   executable="python -m pat.hacc.cinema", 
						   arguments=[ "--input-file", arg1 , "--output-file" , arg2 ],
						   configurations=configurations,
						   environment=environment)

		cinema_job.add_command("cd " + self.json_data["foresight-home"] + "/Analysis/")

		#cinema_job = j.Job(name="cinema",
		#				   execute_dir="cinema",
		#				   executable=cinema_exe,
		#				   arguments=["--input-file", self.workflow_dir + "/" + self.json_path],
		#				   configurations=configurations,
		#				   environment=environment)


		# make dependent on all previous workflow jobs and add to workflow
		#cinema_job.add_parents(*self.jobs)
		#self.add_job(cinema_job)

		if ( len(self.jobs) > 0):
			#print(self.jobs)
			for job in self.jobs:
				cinema_job.add_parents(job)
			self.add_job(cinema_job)
		else:
			self.add_job(cinema_job)



def main():
	# parse command line
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("--input-file", required=True)
	parser.add_argument("--analysis-cinema", action="store_true", help="run analysis and cinema job only")
	parser.add_argument("--cbench",          action="store_true", help="run cbench only")
	parser.add_argument("--cinema",          action="store_true", help="run cinema only")
	parser.add_argument("--preview", 		 action="store_true", help="preview the job, create scripts, ... but won't run")
	opts = parser.parse_args()


	# read input JSON file
	wflow_data = futils.read_json(opts.input_file)


	# create Workflow instance
	wflow_dir = wflow_data["project-home"] + wflow_data["wflow-path"]
	wflow = HACCWorkflow("wflow", wflow_data, workflow_dir=wflow_dir)


	# make directories
	if not os.path.exists(wflow_dir + "/cbench"):
		os.makedirs(wflow_dir + "/cbench")
	if not os.path.exists(wflow_dir + "/halo"):
		os.makedirs(wflow_dir + "/halo")
	if not os.path.exists(wflow_dir + "/spectrum"):
		os.makedirs(wflow_dir + "/spectrum")
	if not os.path.exists(wflow_dir + "/cinema"):
		os.makedirs(wflow_dir + "/cinema")


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
