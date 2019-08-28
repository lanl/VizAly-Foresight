#! /usr/bin/env python

import sys 
import argparse
import os
import csv
import operator

from pat.utils import file_utilities as futils
from pat.utils import plot_utilities as putils
from pat.utils import cinema
from pat.utils import job as j




class nyx_cinema(cinema.CinemaWorkflow):

	def prepare_cinema(self):
		# Open CSV file
		if "metrics-csv" in self.json_data["input"]:
			metrics_csv = self.json_data["input"]["metrics-csv"]
		else:
			metrics_csv = self.json_data["project-home"] +  self.json_data['wflow-path'] + "/cbench/" + self.json_data['cbench']['output']['metrics-file'] + ".csv"

		output_file_name = self.json_data["project-home"] +  self.json_data['wflow-path'] + "/cinema/" + "data.csv"


		all = []
		#reader = futils.open_csv_file(metrics_csv)
		with open(metrics_csv,'r') as csvinput:
			reader = csv.reader(csvinput)

			# Modify Cinema files
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

			
			futils.write_csv(output_file_name, all)


	# Converts a string to an operator
	def validate(self, operand, relate, result):
		ops = { "<": operator.lt }
		return ops[relate]( abs(operand-1.0), result)

	# Process checks found in "cinema-plots" "plotting" "checks"
	def is_valid(self, pk_ratio):
		valid = True
		for check in self.json_data['cinema-plots']['plotting']['checks']:
			for item in pk_ratio:
				if not self.validate(item, check['operator'],  check['result']):
					valid = False 
					break
			
		return valid

				

	def create_plots(self):
		output_path = self.json_data['project-home'] + self.json_data['wflow-path']
		output_plot_path = output_path + "/plots"
		x_range = self.json_data['cinema-plots']['plotting']['x-range']


		for ana in self.json_data['pat']['analysis']:

			plot_title = ana['title']
			to_plot = []  # all the items to plot

			k_list = []
			orig_pk = []

			# Find the original file
			for file in ana['files']:
				if (file['name']=="orig"):
					k_list  = futils.extract_csv_col(file['path'], ' ', 2)
					orig_pk = futils.extract_csv_col(file['path'], ' ', 3)

			# Process the other files
			for file in ana['files']:
				if (file['name']!="orig"):
					temp_pk = futils.extract_csv_col(file['path'], ' ', 3)
					if (temp_pk is not None):
						pk_ratio = [i / j for i, j in zip(temp_pk, orig_pk)]
						this_tuple = (pk_ratio, file['name']) #array, name

						# Check if passes test
						if "checks" in self.json_data["cinema-plots"]["plotting"]:
							if self.is_valid(pk_ratio):
								to_plot.append(this_tuple)
						else:
							to_plot.append(this_tuple)

			putils.plotScatterGraph(k_list, 'k', 'pk-ratio', plot_title, output_plot_path, x_range, to_plot)



# Parse Input
parser = argparse.ArgumentParser()
parser.add_argument("--input-file")
opts = parser.parse_args()


# Create Cinema DB
cinema = nyx_cinema( opts.input_file )
cinema.create_plots()
cinema.create_cinema()
