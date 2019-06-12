#! /usr/bin/env python

import sys, os, json, csv


def list_files_in_folder(folder):
	files = []

	# r=root, d=directories, f = files
	for r, d, f in os.walk(folder):
		for file in f:
			files.append(file)

	return files



def create_folder(path):
	try:  
		os.mkdir(path)
	except OSError:  
		print ("Creation of the directory %s failed" % path)
	else:  
		print ("Successfully created the directory %s " % path)


def validate_num_args(num_arguments, min_count):
	if num_arguments < min_count:
		print ("At least " + str(min_count) + " arguments are needed!")
		exit(1)

	return True


def splitString(filename, char):
	k = filename.rfind(char)
	return filename[:k+1], filename[k+1:]


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

		except ValueError, e:
			print ("CSV file " + filename + " is invalid! " + e)
			exit(1)
	
		


def read_json(filename):
	# Check if file exists
	if ( not os.path.isfile(filename) ):
		print( filename + " does not exist! ")
		return None
		
	# Open json file
	with open(filename, "r") as read_file:
		try:
			json_data = json.load(read_file)
			return json_data
		except ValueError, e:
			print ("Json file " + filename + " is invalid! " + e)
			return None


