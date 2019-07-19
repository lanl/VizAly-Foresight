#!/usr/bin/env python
"""
To run executable do:
python Analysis/selective_analysis.py --input-file inputs/analysis/selective_analysis.json
"""

import argparse, os, csv

from pat import file_utilities as futils
from pat import plot_utilities as putils
from pat import cinema
from pat import Job as j


class selective_analysis(cinema.CinemaWorkflow):
	"""Allows the creatrion of specific cinema databases based on user input"""

	def __init__(self, json_data):
		self.json_data = json_data

		# Specify names / paths
		self.output_path = self.json_data["project-home"] + "/" + self.json_data["wflow-path"]
		self.output_metrics_path = self.output_path + "/cbench/"
		self.output_plot_path    = self.output_path + "/plots/"

		futils.create_folder(self.output_path)	# Create folder for output path
		

	def create_metrics(self):
		futils.create_folder(self.output_metrics_path) # Create the folder to store the CSV file
		self.cinema_metrics_csv = self.output_metrics_path + "/metrics_.csv"

		rows = []
		csv_header = self.json_data["pat"]["metrics-header"]
		rows.append( [csv_header] )


		# Loop through analysis entries in inout file
		for ana in self.json_data["pat"]["analysis"]:
			title = ana["title"]

			# Fill in data from each metric file
			for entry in ana["files"]:
				if (entry["name"] != "orig"):
					input_csv_file = entry["csv-file"]

					with open(input_csv_file,'r') as csv_input:
						reader = csv.reader(csv_input)

						for csv_row in reader:
							for field in entry["compressor_field_param"]:
								if csv_row[0] == field:
									csv_row.insert(1,title)	# Specify the parameter this file is about
									rows.append(csv_row)	# add the rest of the rows

		# write the data
		with open(self.cinema_metrics_csv, 'w') as csvoutput:
			writer = csv.writer(csvoutput, escapechar=' ', quoting=csv.QUOTE_NONE)
			writer.writerows(rows)


	def create_plots(self):
		futils.create_folder(self.output_plot_path) # Create the plots folder
		x_range = self.json_data["cinema-plots"]["plotting"]["x-range"]

		for ana in self.json_data["pat"]["analysis"]:
			plot_title = ana["title"]

			to_plot = []  # all the items to plot

			k_list = []
			orig_pk = []

			# Find the original file
			for entry in ana["files"]:
				if (entry['name'] == "orig"):
					k_list  = futils.extract_csv_col(entry["pwr-spectrum"], ' ', 2)
					orig_pk = futils.extract_csv_col(entry["pwr-spectrum"], ' ', 3)

			for entry in ana["files"]:
				if (entry['name'] != "orig"):

					temp_pk = futils.extract_csv_col(entry["pwr-spectrum"], ' ', 3)
					if (temp_pk is not None):
						pk_ratio = [i / j for i, j in zip(temp_pk, orig_pk)]
						this_tuple = (pk_ratio, entry["name"]) #array, name
						to_plot.append(this_tuple)

			putils.plotScatterGraph(k_list, 'k', 'pk', plot_title, self.output_plot_path, x_range, to_plot)



	def prepare_cinema(self):
		metrics_csv 	 = self.cinema_metrics_csv
		output_file_name = self.output_metrics_path + "/data.csv"

		all = []
		with open(metrics_csv,'r') as csvinput:
			reader = csv.reader(csvinput)

			# Modify Cinema files
			row = next(reader)
			row.append('FILE')
			all.append(row)

			for row in reader:
				row.append( row[1] + ".png")
				all.append(row)

			futils.write_csv(output_file_name, all)



# Parse Input
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--input-file")
opts = parser.parse_args()

print (opts.input_file)

# Read input JSON file
wflow_data = futils.read_json(opts.input_file)

sa = selective_analysis(wflow_data)
sa.create_metrics()
sa.create_plots()
sa.create_cinema()