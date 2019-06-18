#! /usr/bin/env python
import sys, os, json, csv, argparse
from collections import OrderedDict


# Utils
def validate_num_args(num_arguments, min_count):
	if num_arguments < min_count:
		print ("At least " + str(min_count) + " arguments are needed!")
		exit(1)

	return True


def splitString(filename, char):
	k = filename.rfind(char)
	return filename[:k+1], filename[k+1:]



# Commands
def run_command(cmd):
	return subprocess.call(cmd, shell=True)


# Folders
def create_folder(path):
	try:  
		os.mkdir(path)
	except OSError:  
		print ("Creation of the directory %s failed" % path)
	else:  
		print ("Successfully created the directory %s " % path)


def list_files_in_folder(folder):
	files = []

	# r=root, d=directories, f = files
	for r, d, f in os.walk(folder):
		for file in f:
			files.append(file)

	return files



# CSV
def open_csv_file(csv_file):
	with open(csv_file,'r') as csvinput:
		reader = csv.reader(csvinput)
		return reader


def write_csv(csv_file, contents):
	with open(csv_file, 'w') as csvoutput:
		writer = csv.writer(csvoutput, lineterminator='\n')
		writer.writerows(contents)


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
			print ("CSV file " + filename + " is invalid! " + e)
			exit(1)


# JSON
def write_json_file(json_file_name, contents):
	with open(json_file_name, 'w') as file_out:
		file_out.write( json.dumps(contents, indent=4) )
		file_out.close()


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
			print ("Json file " + filename + " is invalid! " + e)
			return None
