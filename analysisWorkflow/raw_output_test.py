#!/usr/bin/python

import sys, json, os
from gioSqlite import *


def output_file(file, directory, result):
	if not os.path.exists(directory):
		os.makedirs(directory)

	run_script_file = open((directory + "/" + file), "w") 
	run_script_file.write(result)
	run_script_file.close()



def analysis(query_mgr):
	intput_files = [ 
		"/projects/groups/vizproject/dssdata/cosmo/Argonne_L360_HACC001/STEP499/m000.full.mpicosmo.499", 
		"/projects/exasky/data/__orig__m000.full.mpicosmo.499"
	]

	count = 0
	for entry in intput_files:
		print("Running analysis on: " + entry)
		
		# Load dataset
		table_name = "foo_" + str(count)
		query_mgr.createTable(table_name, entry )

		# Query
		query = "select count(*) from __TABLE__"
		query = query.replace("__TABLE__", table_name)

		# Run query
		result = query_mgr.runQueryOutputString(query)

		# Save output
		filename = (entry.split("/")[-1])	# extract the filename
		result_filename = filename + '_' + str(count)

		output_file(result_filename, "raw_results", result)

		count = count + 1



if __name__ == "__main__":
	gio_sqlite_lib = "/projects/exasky/VizAly_genericio/genericio/frontend/GenericIOSQLite.so"

	# Initialize query tool
	query_mgr = GioSqlite3()
	if query_mgr.loadGIOSqlite(gio_sqlite_lib) == -1:
		exit()


	# Run analysis
	analysis(query_mgr)


# Run as:
# python workflow_b.py workflow_input.json

# Results will be in subdirectory results