#!/usr/bin/python

import sys, os, json, csv
from collections import OrderedDict

import file_utilities as f
import plot_utilities as p
import job as j


class Analysis_Workflow:
	sim = ""
	original_json_filename = ""
	analysis_json_filename = ""


	def __init__(self, sim_name):
		self.sim = sim_name


	# Append input generated from CBench to input file
	def init_analysis(self, num_input_params, json_file):
		# Check number of arguments
		self.json_filename = json_file
		status = f.validate_num_args(num_input_params, 2)
		if not status:
			return False

		# Check Validity
		self.original_json_data = f.read_json(self.json_filename)
		if self.original_json_data == None:
			return False


		# Open File to append "input-files"
		input_json_filename = self.json_filename
		output_json_filename = self.original_json_data['simulation-analysis']['analysis-folder'] + ".json"

		with open(input_json_filename, 'r') as file_in:
			json_data = json.load(file_in, object_pairs_hook=OrderedDict)


		# Create input settings
		orig_path_filename = f.splitString(json_data['input']['filename'],'/')

		orig_item = {
			"output-prefix" : "orig",
			"path" : json_data['input']['filename']
		}
		json_data['simulation-analysis']['input-files'].append(orig_item)


		for _file in json_data['compressors']:
			json_item = {
				"output-prefix" : _file["output-prefix"],
				"path" : json_data['output']['run-path'] +  json_data['output']['output-decompressed-location'] + "/" + _file['output-prefix'] + "__" + orig_path_filename[1]
			}

			json_data['simulation-analysis']['input-files'].append(json_item)


		# Write out the file
		with open(output_json_filename, 'w') as file_out:
			file_out.write( json.dumps(json_data, indent=4) )
			file_out.close()

		self.analysis_json_filename = output_json_filename
		self.json_data = f.read_json(self.analysis_json_filename)


	def run_analysis(self):
		pass

	def create_plots(self):
		pass

	def create_CinemaDB(self):
		pass