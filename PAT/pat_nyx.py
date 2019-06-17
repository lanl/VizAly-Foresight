#! /usr/bin/env python

import argparse
import os
from pat import file_utilities as futils
from pat import plot_utilities as putils
from pat import workflow


def write_analysis_input(analysis_json_filename):
	temp_json_data =  futils.read_json(analysis_json_filename)
	
	# Create analysis settings
	currentPath = os.getcwd()

	for ana in temp_json_data['simulation-analysis']['analysis-tool']['analytics']:
		for item in ana['type']:
			json_item = {
				"title" : ana['name'] + "_" + item,
				"files" : []
			}

			for inputItem in temp_json_data['simulation-analysis']['input-files']:
				input_item = { 
					'name' : inputItem["output-prefix"],
					'path' : currentPath + "/" + temp_json_data['simulation-analysis']['analysis-folder'] + "/" + ana['name'] + "/" + inputItem['output-prefix'] + item + ana['postfix']
				}

				json_item['files'].append(input_item)

			
			temp_json_data['simulation-analysis']['analysis'].append(json_item)

	futils.write_json_file(analysis_json_filename, temp_json_data)

	return analysis_json_filename
	

def run_analysis(json_file):
	json_data = futils.read_json(json_file)

	# Set some environment variables
	os.system("cd " + json_data['simulation-analysis']['analysis-tool']['gimlet-home'])


	# Create Folder for analytics
	path = os.getcwd()
	path = path + "/" + self.json_data['simulation-analysis']['analysis-folder']
	if not os.path.exists(path):
		os.mkdir(path)


	# Run analytics
	for sim_analytics in self.json_data['simulation-analysis']['analysis-tool']['analytics']:
		for file in self.json_data['simulation-analysis']['input-files']:
			# Run gimlet analysis
			analysis_path = path + "/" + sim_analytics['name']
			if not os.path.exists(analysis_path):
				os.mkdir(analysis_path)

			cmd = "mpirun -np 16 " + sim_analytics['path'] + " " + file['path'] + " " + analysis_path + "/" + file['output-prefix']
			os.system(cmd)


def create_cinema(json_file):
	# Read Json file
	json_data = futils.read_json(json_file)

	# Open CSV file
	metrics_csv = self.json_data['output']['run-path'] + self.json_data['output']['metricsfname'] + ".csv "
	reader = futils.open_csv_file(metrics_csv)

	# Modify Cinema files
	all = []
	row = next(reader)
	row.append('FILE_SimStats_Pk')
	row.append('FILE_lya_all_axes_x_Pk')
	row.append('FILE_lya_all_axes_y_Pk')
	row.append('FILE_lya_all_axes_z_Pk')
	all.append(row)


	values = ["sim_stats_rhob.png", "sim_stats_rhodm.png", "sim_stats_temp.png", "sim_stats_velmag.png", "sim_stats_velmag.png", "sim_stats_vz.png"]
	count = 0
	for row in reader:
		row.append(values[count])
		row.append("lya_all_axes_x.png")
		row.append("lya_all_axes_y.png")
		row.append("lya_all_axes_z.png")
		all.append(row)

		count = count + 1
		if (count == 6):
			count = 0

	futils.write_csv("data.csv", all)


	# Create cinema database
	cinema_database = self.json_data['simulation-analysis']['analysis-folder'] + ".cdb"
	img_files = list_files_in_folder( self.json_data['simulation-analysis']['analysis-folder'] )
	self.create_CinemaDB(cinema_database, "data.csv", img_files)


def create_plots(self):
	# Read Json file
	json_data = futils.read_json(json_file)

	path = self.json_data['simulation-analysis']['analysis-folder']
	csv_file_path = self.json_data['output']['run-path'] + self.json_data['output']['metricsfname'] + ".csv"
	x_range = self.json_data['simulation-analysis']['plotting']['x-range']

	for ana in self.json_data['simulation-analysis']['analysis']:
		plot_title = ana['title']
		to_plot = []  # all the items to plot

		k_list = []
		orig_pk = []

		# Find the original file
		for file in ana['files']:
			if (file['name']=="orig"):
				k_list  = f.extract_csv_col(file['path'], ' ', 2)
				orig_pk = f.extract_csv_col(file['path'], ' ', 3)

		for file in ana['files']:
			if (file['name']!="orig"):
				print (file['path'])

				temp_pk = f.extract_csv_col(file['path'], ' ', 3)
				if (temp_pk is not None):
					pk_ratio = [i / j for i, j in zip(temp_pk, orig_pk)]
					this_tuple = (pk_ratio, file['name']) #array, name
					to_plot.append(this_tuple)

		putils.plotScatterGraph(k_list, 'k', 'pk', plot_title, path, x_range, to_plot)




# Parse Input
parser = argparse.ArgumentParser()
parser.add_argument("--input-file")
parser.add_argument("--submit", action="store_true")
opts = parser.parse_args()

# Read input JSON file
wflow_data = futils.read_json(opts.input_file)

# create Workflow instance
wflow_dir = wflow_data["project-home"]
wflow = workflow.Workflow("wflow", wflow_data, workflow_dir=wflow_dir)


# add CBench job
wflow.add_cbench_job()
wflow.add_analysis_jobs()
wflow.add_plotting_jobs()


# write submission script
wflow.write_submit()


# submit workflow
if opts.submit:
    wflow.submit()


"""
python pat_nyx.py --input-file ../inputs/nyx/NYX_wflow.json 
"""