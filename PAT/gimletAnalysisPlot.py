#!/usr/bin/python

import sys, json, os, csv
import pandas as pd
from collections import OrderedDict

from analysis_lib import file_utilities as f
from analysis_lib import plot_utilities as p

from analysis_lib.analysis_workflow import Analysis_Workflow



def write_analysis_input(input_json_filename, output_json_filename):
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


	# Create analysis settings
	currentPath = os.getcwd()

	for ana in json_data['simulation-analysis']['analysis-tool']['analytics']:
		for item in ana['type']:
			json_item = {
				"title" : ana['name'] + "_" + item,
				"files" : []
			}

			for inputItem in json_data['simulation-analysis']['input-files']:
				input_item = { 
					'name' : inputItem["output-prefix"],
					'path' : currentPath + "/" + json_data['simulation-analysis']['analysis-folder'] + "/" + ana['name'] + "/" + inputItem['output-prefix'] + item + ana['postfix']
				}

				json_item['files'].append(input_item)

			
			json_data['simulation-analysis']['analysis'].append(json_item)


	with open(output_json_filename, 'w') as file_out:
		file_out.write( json.dumps(json_data, indent=4) )
		file_out.close()


def createCommand(hdfpath, input_file, output_file, field):
	command = hdfpath + "h5copy" + " -v -i \"" + input_file + "\" -o \"" + output_file + "\" -s \"" + field + "\" -d \"" + field + "\""
	return command


def runGimlet(json_data):
	# Set some environment variables
	os.system("source " + json_data['simulation-analysis']['evn_path'])
	os.system("cd " + json_data['simulation-analysis']['analysis-tool']['gimlet-home'])

	# Find original file
	index_orig = 0
	for file in json_data['simulation-analysis']['input-files']:
		if (file["output-prefix"] != "orig"):
			break
		
	# Copy uncompressed attributes from original to compressed; needed for gimlet
	for file in json_data['simulation-analysis']['input-files']:
		if (file["output-prefix"] != "orig"):
			# Copy missing attributes
			hdf5_path = json_data['simulation-analysis']['analysis-tool']['hdf5_path']

			# field group
			command = createCommand(hdf5_path, json_data['simulation-analysis']["input-files"][index_orig]["path"], file["path"], "/domain")
			os.system(command)

			# universe group
			command = createCommand(hdf5_path, json_data['simulation-analysis']["input-files"][index_orig]["path"], file["path"], "/universe")
			os.system(command)


	# Create Folder for analytics
	path = os.getcwd()
	path = path + "/" + json_data['simulation-analysis']['analysis-folder']
	if not os.path.exists(path):
		os.mkdir(path)


	# Run analytics
	for sim_analytics in json_data['simulation-analysis']['analysis-tool']['analytics']:
		for file in json_data['simulation-analysis']['input-files']:
			# Run gimlet analysis
			analysis_path = path + "/" + sim_analytics['name']
			if not os.path.exists(analysis_path):
				os.mkdir(analysis_path)

			cmd = "mpirun -np 16 " + sim_analytics['path'] + " " + file['path'] + " " + analysis_path + "/" + file['output-prefix']
			os.system(cmd)


def gimletPlot(json_data):
	path = json_data['simulation-analysis']['analysis-folder']

	for ana in json_data['simulation-analysis']['analysis']:
		plot_title = ana['title']

		to_plot = []  # all the items to plot

		k_list = []
		orig_pk = []
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

		p.plotGraph(k_list, 'k', 'pk', plot_title, path, to_plot)



def createCinemaDB(db_name, json_data):
	# Create cinema db folder
	cinema_database = json_data['simulation-analysis']['analysis-folder'] + ".cdb"
	print cinema_database
	f.create_folder(cinema_database)

	# Copy files to cinema
	cmd = "cp " + json_data['simulation-analysis']['analysis-folder'] + "/*.png " + cinema_database
	os.system(cmd)


	# Copy metrics to cinema
	cmd = "cp " + json_data['output']['run-path'] + json_data['output']['metricsfname'] + ".csv " + cinema_database + "/data.csv"
	os.system(cmd)


	# Modify csv to create cinema database
	cinema_csv = cinema_database + "/data.csv"


	with open(cinema_csv,'r') as csvinput:
		reader = csv.reader(csvinput)

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


	with open(cinema_csv, 'w') as csvoutput:
		writer = csv.writer(csvoutput, lineterminator='\n')
		writer.writerows(all)
	


if __name__ == "__main__":
	analysis = Analysis_Workflow("NYX")
	analysis.init_analysis( len(sys.argv), sys.argv[1] )

	# create analysis input file
	analysis_filename = analysis.json_data['simulation-analysis']['analysis-folder'] + ".json"
	write_analysis_input(sys.argv[1], analysis_filename)
	analysis.set_analysis(analysis_filename)


	# Run the analysis
	#runGimlet(analysis.json_data)

	# Draw Graph
	gimletPlot(analysis.json_data)
	
	# Create Cinema DB
	createCinemaDB("test", analysis.json_data)

"""
python gimletAnalysisPlot.py ../inputs/nyx_all.json
"""