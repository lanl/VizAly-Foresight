#! /usr/bin/env python

import sys 
import argparse
import os
import csv
import operator

from draw.utilities import utilities as utils
from draw.plot_utilities import plot_utilities as putils


def create_plots(output_plot_path, ):
	#output_path = self.json_data['project-home'] + self.json_data['wflow-path']
	#output_plot_path = output_path + "/plots"

	#has_range = False;
	#if "x-range" in self.json_data['cinema-plots']['plotting']:
	#	x_range = self.json_data['cinema-plots']['plotting']['x-range']
	#	has_range = True



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


		# Check range limit
		range_count = 0
		for x in k_list:
			#if has_range == True:
			#	if (x > x_range[1]):
			#		break
			#else:
			#	has_range = 1000000

			range_count =  range_count + 1

		# Process the other files
		for file in ana['files']:
			if (file['name']!="orig"):
				temp_pk = futils.extract_csv_col(file['path'], ' ', 3)
				if (temp_pk is not None):
					pk_ratio = [i / j for i, j in zip(temp_pk, orig_pk)]
					this_tuple = (pk_ratio, file['name']) #array, name

					# Check if passes test
					if "checks" in self.json_data["cinema-plots"]["plotting"]:
						if self.is_valid(pk_ratio, range_count):
							to_plot.append(this_tuple)
					else:
						to_plot.append(this_tuple)

		putils.plotScatterGraph(k_list, 'k', 'pk-ratio', plot_title, output_plot_path, x_range, to_plot)