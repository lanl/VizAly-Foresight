#!/usr/bin/python

import csv
import os
import sys
import json
import matplotlib
import itertools
from collections import OrderedDict

matplotlib.use('agg')
import matplotlib.pyplot as plt

def read_json(filename):
	# Check if file exists
	if ( not os.path.isfile(filename) ):
		print( filename + " does not exist! ")
		return None
		
	# Open json file
	with open(filename, "r") as read_file:
		try:
			json_data = json.load(read_file,  object_pairs_hook=OrderedDict)
			return json_data
		except ValueError as e:
			print ("JSON file " + filename + " is invalid! \n\t" + str(e) + "\n")
			return None


# From: https://stackoverflow.com/questions/29461608/matplotlib-fixing-x-axis-scale-and-autoscale-y-axis
def autoscale_y(ax,margin=0.1):
	import numpy as np

	def get_bottom_top(line):
		xd = line.get_xdata()
		yd = line.get_ydata()
		lo,hi = ax.get_xlim()
		y_displayed = yd[((xd>lo) & (xd<hi))]
		h = np.max(y_displayed) - np.min(y_displayed)
		bot = np.min(y_displayed)-margin*h
		top = np.max(y_displayed)+margin*h
		return bot,top

	lines = ax.get_lines()
	bot,top = np.inf, -np.inf

	for line in lines:
		new_bot, new_top = get_bottom_top(line)
		if new_bot < bot: bot = new_bot
		if new_top > top: top = new_top

	ax.set_ylim(bot,top)


def plotScatterGraph(x, x_label, y_label, title, path, x_range, list_of_tuples):
	fig = plt.figure()
	ax = plt.subplot(1,1,1)

	plt.title(title)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.grid(True)

	marker = itertools.cycle(('o', 's', '*', 'P', 'D', 'X', 'v', 'H', '2', '+'))

	for item in list_of_tuples:
		ax.semilogx(x, item[0], label=item[1], marker=next(marker))
		ax.legend(bbox_to_anchor=[1.05,1],loc='upper left',borderaxespad=0)

	plt.xlim(x_range[0], x_range[1])

	fig = plt.gcf()
	autoscale_y(ax)
	plt.show()
	fig.savefig(path +'/' + title + '.png', dpi=100, bbox_inches="tight")



def extract_csv_col(filename, delimiter_char, colpos):
	if ( not os.path.isfile(filename) ):
		print( filename + " does not exist! ")
		return None

	with open(filename) as csv_file:
		try:
			csv_data = csv.reader(csv_file, delimiter=delimiter_char)
			col = []
			for row in csv_data:
				col.append( float(row[colpos]) )

			return col

		except ValueError as e:
			print ("CSV file " + filename + " is invalid! " + str(e))
			exit(1)


def create_plots(filename):
	json_data = read_json(filename)
 
	output_path = json_data['project-home'] + json_data['wflow-path']
	output_plot_path = output_path + "/plots"


	# has_range = False;
	# if "x-range" in json_data['cinema-plots']['plotting']:
	# 	x_range = json_data['cinema-plots']['plotting']['x-range']
	# 	has_range = True
	has_range = True
	x_range = [0,11]



	for ana in json_data['analysis']['output-files']:

		plot_title = ana['title']
		to_plot = []  # all the items to plot

		k_list = []
		orig_pk = []

		# Find the original file
		for file in ana['files']:
			if (file['name']=="orig"):
				print(file['name'])
				k_list  = extract_csv_col(file['path'], ' ', 2)
				orig_pk = extract_csv_col(file['path'], ' ', 3)
				print(k_list)
				print(orig_pk)


		# Check range limit
		range_count = 0
		for x in k_list:
			if has_range == True:
				if (x > x_range[1]):
					break
			else:
				has_range = 1000000

			range_count =  range_count + 1

		# Process the other files
		for file in ana['files']:
			if (file['name']!="orig"):
				temp_pk = extract_csv_col(file['path'], ' ', 3)
				if (temp_pk is not None):
					pk_ratio = [i / j for i, j in zip(temp_pk, orig_pk)]
					this_tuple = (pk_ratio, file['name']) #array, name

					# Check if passes test
					#if "checks" in json_data["cinema-plots"]["plotting"]:
					#	if is_valid(pk_ratio, range_count):
					#		to_plot.append(this_tuple)
					#else:
					#	to_plot.append(this_tuple)
					to_plot.append(this_tuple)

		plotScatterGraph(k_list, 'k', 'pk-ratio', plot_title, output_plot_path, x_range, to_plot)
 

# Create Cinema DB
create_plots( sys.argv[1] )

