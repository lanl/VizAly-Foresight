import sys
import os
import json
import csv
import argparse
from collections import OrderedDict

from pat.utils import utilities as utils
from pat.utils import plot_utilities as putils
from pat.utils import job as j



class Workflow(object):
	""" Class that describes workflow for Slurm scheduler.
	Note that jobs ordering appended to workflow has meaning in this simple API.
	"""

	def __init__(self, name, json_data, workflow_dir=""):
		# store meta-data about the workflow
		self.name = name
		self.json_data = json_data

		self.json_path = self.json_data['project-home'] + self.json_data['wflow-path'] + "/" + self.name + ".json"

		# store a list of jobs in workflow
		self.jobs = []

		# location to store workflow files
		self.workflow_dir = workflow_dir

		# attributes for workflow construction
		self.submit_file = None

		# Add git tag
		git_tag = utils.get_git_version(self.json_data['foresight-home'])
		git_tag = git_tag.strip('\n')
		self.json_data["git-tag"] = git_tag


	def fill_reduc_output_presets(self):
		""" Fill in  ["data-reduction"]["output-files"] with some user specified presets
		"""
		for item in self.json_data["data-reduction"]["pre-compressed"]:
			json_item = {
				"output-prefix" : item["output-prefix"],
				"path" : item["path"]
			}

			self.json_data['data-reduction']['output-files'].append(json_item)


	def fill_reduc_csv(self, filename):
		""" Append cbench csv with other results
		"""
		header = utils.get_csv_header(filename)

		for item in self.json_data["data-reduction"]["pre-compressed"]:
			mylist = [''] * len(header)
   
			for params in item["compression-params"]:
				index = 0
				for h in header:
					if h == params.key():
						mylist[index] = params.value()
					index = index + 1
    
			utils.append_list_as_row(filename, mylist)


	def fill_cbench_output_files(self):
		""" Fill in  ["data-reduction"]["output-files"]
		"""
		base_path = self.json_data['project-home'] + self.json_data['wflow-path']

		## Remove all entries if any - useful for rerunning an existing workflow
		self.json_data['data-reduction']['output-files'].clear()

		# Add the original
		orig_path_filename = utils.splitString(self.json_data['input']['filename'],'/')
		orig_item = {
			"output-prefix" : "orig",
			"path" : self.json_data['input']['filename']
		}
		self.json_data['data-reduction']['output-files'].append(orig_item)

		# Add decompressed ones
		for _file in self.json_data['data-reduction']['cbench-compressors']:
			json_item = {
				"output-prefix" : _file["output-prefix"],
				"path" : base_path + "/reduction/cbench/" + self.json_data['data-reduction']['cbench-output']['output-decompressed-location'] + "/" + _file['output-prefix'] + "__" + orig_path_filename[1]
			}

			self.json_data['data-reduction']['output-files'].append(json_item)


	def add_cbench_job(self):
		""" Adds a CBench job to the workflow.
		"""
	
		# Set up environment
		execute_dir = "reduction/cbench"

		if "configuration" in self.json_data["data-reduction"].keys():
			configurations = list(sum(self.json_data["data-reduction"]["configuration"].items(), ()))
		else:
			configurations = None

		if "evn_path" in self.json_data["data-reduction"]:
			environment =  self.json_data["foresight-home"] + self.json_data["data-reduction"]["evn_path"]
		else:
			environment = None

		# add a single CBench job to workflow for entire sweep
		cbench_job = j.Job(name="cbench",
						 job_type = "reduction",
						 execute_dir=execute_dir,
						 executable=self.json_data["data-reduction"]["path"],
						 arguments=[self.json_path],
						 configurations=configurations,
						 environment=environment)
		cbench_job.add_command("mkdir -p logs")
		self.add_job(cbench_job)
  
		# Fill output files for cbench
		self.fill_cbench_output_files()


	def add_job(self, this_job, dependencies="all_previous", filter=""):
		""" Adds a job to the workflow.
		"""
		if dependencies == "all_previous":
			if len(self.jobs) > 0:
				job_list = []
				for job in self.jobs:
					job_list.append(job)

				this_job.add_parents(job_list)
			self.jobs.append(this_job)


		elif dependencies == "type":
			if len(self.jobs) > 0:
				job_list = []
				for job in self.jobs:
					if job.job_type == filter:
						job_list.append(job)

				this_job.add_parents(job_list)

			self.jobs.append(this_job)


		elif dependencies == "single":
			if len(self.jobs) > 0:
				found_job = False
				for job in self.jobs:
					if job.name == filter:
						this_job.add_parents([job])
						found_job = True
						break

			self.jobs.append(this_job)


		else:
			self.jobs.append(this_job)



	def add_data_reduction_jobs(self):
		""" Adds analysis jobs to workflow that do not produce final products.
		"""
		raise NotImplementedError("Implement the `add_data_reduction` function to your workflow!")


	def add_analysis_jobs(self):
		""" Adds analysis jobs to workflow that do not produce final products.
		"""
		raise NotImplementedError("Implement the `add_analysis` function to your workflow!")


	def add_cinema_plotting_jobs(self):
		""" Adds plotting jobs to workflow that produce final products.
		"""
		raise NotImplementedError("Implement the `add_plotting_jobs` function to your workflow!")




	def write_submit(self):
		""" Writes Slurm workflow files.
		"""

		# get foresight dir
		foresight_home = self.json_data["foresight-home"]

		# create workflow output dir
		# change to workflow output dir
		utils.create_folder(self.workflow_dir)
		os.chdir(self.workflow_dir)
	
		# write JSON data
		if not os.path.exists(os.path.dirname(self.json_path)):
			os.makedirs(os.path.dirname(self.json_path))
		with open(self.json_path, 'w') as fp:
			fp.write( json.dumps(self.json_data, indent=4) )
			fp.close()

		# create submission script
		self.submit_path = self.name + ".sh"
		with open(self.submit_path, "w") as fp:
			fp.write("#! /bin/bash\n")

		# loop over each job
		for i, job in enumerate(self.jobs):
	
			# get a unique name and index
			job.name = job.name if job.name != None else "job_{}".format(i)
			job._idx = i
	
			# figure out dependencies job indices
			parent_idxs = [j._idx for j in job._parents]
			if parent_idxs == []:
				depends_str = ""
			else:
				depends_str = "--dependency=afterok:" + ":".join(["$jid{}".format(j) for j in parent_idxs])
	
			# create workflow directory
			if not os.path.exists(job.execute_dir):
				os.makedirs(job.execute_dir)
	
			# write wrapper script for job
			path = job.execute_dir + "/" + job.name + ".sh"
			with open(path, "w") as fp:
				fp.write("#! /bin/bash\n")
				for key, val in zip(job.configurations[::2], job.configurations[1::2]):
					fp.write("#SBATCH --{}={}\n".format(key, val))

				fp.write("date\n")
				fp.write("mkdir -p {}/{}\n".format(self.workflow_dir, job.execute_dir))
				fp.write("cd {}/{}\n".format(self.workflow_dir, job.execute_dir))

				# Pre-job execution command
				if job.commands != []:
					for _cmd in job.commands:
						fp.write(_cmd + "\n")

				if job.environment != None:
					fp.write("source {}\n".format(job.environment))

				cmd = job.executable + " " + " ".join(map(str, job.arguments)) + "\n"
				cmd = cmd.replace("$foresight-home$", foresight_home)
				fp.write(cmd)

				# post-job execution command
				if job.post_commands != []:
					for _cmd in job.post_commands:
						fp.write(_cmd + "\n")

				fp.write("date\n")

			# append job to controller file
			slurm_out_path = job.execute_dir + "/" + job.name + ".slurm.out"
			with open(self.submit_path, "a") as fp:
				fp.write("\n# {}\n".format(job.name))
				fp.write("jid{}=$(sbatch {} --output {} {})\n".format(job._idx, depends_str, slurm_out_path, path))
				fp.write("jid{}=$(echo $jid{} | rev | cut -f 1 -d ' ' | rev)".format(job._idx, job._idx))
	

	def submit(self):
		""" Submits Slurm workflow.
		"""
		os.system("bash {}".format(self.submit_path))
  
  
  
  # Workflow specific utils
  
def parse_args():
	# Parse Input
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("--input-file")
	parser.add_argument("--preview", 		 	action="store_true", help="preview the job, create scripts, ... but won't run")
	
	# run modes
	#parser.add_argument("--all", action="store_true", help="run full pipeline")

	parser.add_argument("--reduction-analysis", action="store_true", help="run data reduction and analysis only")
	parser.add_argument("--analysis-cinema", 	action="store_true", help="run analysis and cinema job only")

	parser.add_argument("--reduction",  	 	action="store_true", help="run data-reduction only")
	parser.add_argument("--analysis",        	action="store_true", help="run analysis only")
	parser.add_argument("--vis",             	action="store_true", help="run cinema only")
	
	opts = parser.parse_args()
	return opts



def process_input(wflow, opts):
	# add jobs to workflow
	if opts.reduction:
		print("Run data reduction only")
		wflow.add_data_reduction_jobs()
	elif opts.analysis:
		print("Run data analysis only")
		wflow.add_analysis_jobs()
	elif opts.vis:
		print("Run data analysis only")
		wflow.add_cinema_plotting_jobs()

	elif opts.reduction_analysis:
		print("Run data reduction + analysis")
		wflow.add_data_reduction_jobs()
		wflow.add_analysis_jobs()
	elif opts.analysis_cinema:
		print("Run analysis + cinema")
		wflow.add_analysis_jobs()
		wflow.add_cinema_plotting_jobs()

	else:	
		print("Run full: data reduction, analysis, and Cinema")
		wflow.add_data_reduction_jobs()
		wflow.add_analysis_jobs()
		wflow.add_cinema_plotting_jobs()
  
	return wflow