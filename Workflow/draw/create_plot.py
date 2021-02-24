#!/usr/bin/python
import sys
import os
import json
import csv
import subprocess
import fileinput
import argparse

import matplotlib
import itertools
matplotlib.use('agg')

import matplotlib.pyplot as plt

# CSV
def extract_csv_col(filename, delimiter_char, colpos, max_rows=-1):
	if ( not os.path.isfile(filename) ):
		print( filename + " does not exist! ")
		return None

	with open(filename) as csv_file:
		try:
			csv_data = csv.reader(csv_file, delimiter=delimiter_char)
			headers = next(csv_data, None)
			col = []
		   
			for row in csv_data:
				col.append( float(row[colpos]) )
		   

			return col

		except ValueError as e:
			print ("CSV file " + filename + " is invalid! " + str(e))
			exit(1)


def create_img(text, filename):
	fig, ax = plt.subplots() 
   
	fig.text(0.2, 0.4, text, fontsize=30, color="black") 
	ax.set(xlim = (0, 8), ylim = (0, 8)) 

	_filename = filename.replace(' ','_')
	plt.savefig(_filename + ".png")



def generic_plot(x, y_orig, x_label, y_label, name, y_tuple_labels=[], plot_type="linear", axis_type="simple"):
	"""
		axis_type options: linear, x_log, y_log, log
		plot_type:  simple (one y value)
					two (plot both y values)
					diff (multiple y values: y_0 - y_1)
					ratio (multiple y values: y_0/y_1)
	"""

	# Extract Data
	num_files = len(y_orig)


	# Process data if needed
	y = []
	if num_files > 1:
		if plot_type == "linear":
			y = y_orig

		elif plot_type == "two":
			for i in range(num_files):
				y.append( y_orig[i] )

		elif plot_type == "diff":
			
			for i in range(num_files-1):
				temp_y = []
				for j in range( len(y_orig[0]) ):
					temp_y.append(  y_orig[i+1][j] - y_orig[0][j] )
				y.append( temp_y )

				
		elif plot_type == "ratio":
			
			for i in range(num_files-1):
				temp_y = []
				for j in range( len(y_orig[0]) ):
					if ( y_orig[i+1][j] == 0):
						temp_y.append( (y_orig[i+1][j]+1) / (y_orig[0][j]+1) )
					else:
						temp_y.append( y_orig[i+1][j] / y_orig[0][j] )
				y.append( temp_y )
	else:
		y = y_orig


	# Create layout
	color_list = ["blue", "red", "orange", "green",  "black"]

	fig = plt.figure()
	ax = plt.subplot(1,1,1)

	plt.title(name)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.grid(True)

	
	for i in range(len(y)):
		line_label=""
		try:
			line_label = y_tuple_labels[i]
		except:
			line_label=""
		plt.plot(x, y[i], '-', color=color_list[i], label=line_label, marker='')

	if axis_type == "x_log":
		plt.xscale('log')
	elif axis_type == "y_log":
		plt.yscale('log')
	elif axis_type == "log":
		plt.xscale('symlog')
		plt.yscale('symlog')
	
	plt.legend()
	filename = name.replace(' ','_')
	plt.savefig(filename + ".png")



def read_files(filenames, seperator=',', x_col=0, y_col=1, limit_type="None", limit_val=1):
	"""
		limit_type: 
			None: read all
			x_less: x less than limit_value
	"""

	num_files = len(filenames)

	x = extract_csv_col(filenames[0], seperator, x_col)
	
	if x == None:
		x = []
		y = []
		return x, y
	

	y_orig = []
	for i in range(num_files):
		y_orig.append( extract_csv_col(filenames[i], seperator, y_col) )


	if limit_type == "x_less":
		count = 0
		for i in x:
			if i < limit_val:
				count = count + 1

		new_x = x[:count]
		new_y = []
		for i in range(num_files):
			new_y.append(y_orig[i][:count])

		return new_x, new_y
	else:
		return x, y_orig



def create_plot(filenames, x_label, y_label, name, seperator=',', x_col=0, y_col=1, limit_type="None", limit_val=1, y_tuple_labels=[], plot_type="linear", axis_type="simple"):
	x, y = read_files(filenames, seperator, x_col, y_col, limit_type, limit_val)
	if x == [] or y == []:
		print("Could not create", name)      
		create_img("Can't create\n"+name, name)
	else:
		generic_plot(x, y, x_label, y_label, name, y_tuple_labels, plot_type, axis_type)



if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-f", "--filenames", nargs="+",  help="input files")
	parser.add_argument("-x", "--x_label", help="label for x-axis")
	parser.add_argument("-y", "--y_label", help="label for y-axis")
	parser.add_argument("-n", "--name", help="title of the plot")

	parser.add_argument("-s", "--seperator", default=',', help="seperator in input files")
	parser.add_argument("-c", "--x_col", default=0, type=int, help="column to extract x value from")
	parser.add_argument("-d", "--y_col", default=1, type=int, help="column to extract y value from")
	parser.add_argument("-l", "--limit", default="None", help="Stop reading after a condigiton is met; Options are: None, x_less")
	parser.add_argument("-v", "--limit_val", default=1, type=float, help="value for option limit")
	parser.add_argument("-t", "--y_tuple_labels", nargs="+", default=[], help="names of the different plots")
	parser.add_argument("-p", "--plot_type", default="linear", help="plot_type - Options are:  simple (one y value),two (plot both y values), diff (multiple y values: y_0 - y_1), ratio (multiple y values: y_0/y_1)")
	parser.add_argument("-a", "--axis_type", default="simple", help="axis type - Options are: linear, x_log, y_log, log")

	args = parser.parse_args()

	create_plot(args.filenames, args.x_label, args.y_label, args.name, args.seperator, args.x_col, args.y_col, args.limit, args.limit_val, args.y_tuple_labels, args.plot_type, args.axis_type)

# python graphDrawing.py -f "data/pw_spec_sz.txt" -x "k"  -y "p(k)" -n "Power Spectrum SZ" --seperator $'\t' --limit "x_less" --limit_val 10 --axis_type "log"
# python graphDrawing.py -f "data/halo_mass_orig.txt" "data/halo_mass_sz01.txt" -x "Halo_mass"  -y "count" -n "SZ Mass Compared" --plot_type "two"  --y_tuple_labels "orig" "SZ" --axis_type "log"
