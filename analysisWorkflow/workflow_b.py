#!/usr/bin/python

import sys, json, os
from gioSqlite import *


def output_file(file, directory, result):
	if not os.path.exists(directory):
		os.makedirs(directory)

	run_script_file = open((directory + "/" + file), "w") 
	run_script_file.write(result)
	run_script_file.close()


def analysis(query_mgr, json_data):
	output_files = []
	halo_file_path = json_data["halo-output-path"]

	count = 0
	for entry in json_data["analysis"]:

		# Query each file
		for halo_file in entry["halo-files"]:
			count = count + 1
			print("Running analysis on: " + halo_file_path + halo_file)

			# Load dataset
			table_name = "foo_" + str(count)
			query_mgr.createTable(table_name, (halo_file_path + halo_file) )

			# Execute Queries
			query_count = 0
			for sql_query in entry["queries"]:

				filename = (halo_file.split("/")[-1])	# extract the filename
				result_filename = filename + '_' + str(query_count)

				query = sql_query.replace("__TABLE__", table_name)
				result = query_mgr.runQueryOutputString(query)

				# Saving the output
				output_file(result_filename, "results", result)
				output_files.append(result_filename)

				query_count = query_count + 1

	return output_files



if __name__ == "__main__":
	# check if the json file is here
	if len(sys.argv) < 2:
		print("Json file needed")
		exit()

	# parse Json file
	with open(sys.argv[1], "r") as read_file:
		json_data = json.load(read_file)


	# Initialize query tool
	query_mgr = GioSqlite3()
	if query_mgr.loadGIOSqlite(json_data["analysis-prams"]["sqlite-path"]) == -1:
		exit()


	# Run analysis
	output_files = analysis(query_mgr, json_data)


# Run as:
# python workflow_b.py workflow_input.json

# Results will be in subdirectory results