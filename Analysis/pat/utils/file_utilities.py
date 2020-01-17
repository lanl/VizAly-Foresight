#! /usr/bin/env python
import sys
import os
import re
import json
import csv
import argparse
import shutil
import subprocess
from collections import OrderedDict


# Utils
def validate_num_args(num_arguments, min_count):
	if num_arguments < min_count:
		print ("At least " + str(min_count) + " arguments are needed!")
		exit(1)

	return True


def parseList(s):
	# s should be in the format ['aa','bb','cc']
	try:
		s1 = s.replace('[','')
		s2 = s1.replace(']','')
		s3 = s2.replace('\'','')

		l = list(map(str, s3.split(',')))
		return l

	except:
		raise argparse.ArgumentTypeError("Something is wrong with the list")


def splitString(filename, char):
	k = filename.rfind(char)
	return filename[:k+1], filename[k+1:]


# From: https://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
def sorted_nicely(l):
	""" Sort the given iterable in the way that humans expect.""" 
 
	convert = lambda text: int(text) if text.isdigit() else text 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(l, key = alphanum_key)



# Commands
def run_command(cmd):
	return subprocess.call(cmd, shell=True)


def get_git_version(path):
	git_tag = os.popen("cd " + path + "; git describe --tags").read()
	return git_tag


# Folders
def create_folder(path):
	if not os.path.exists(path):
		try:  
			os.mkdir(path)
		except OSError:  
			print ("Creation of the directory %s failed" % path)
		else:  
			print ("Successfully created the directory %s " % path)


def delete_folder(path):
	if not os.path.exists(path):
		try:  
			shutil.rmtree(path)
		except OSError:  
			print ("Deletion of the directory %s failed" % path)
		else:  
			print ("Successfully deleted the directory %s " % path)


def list_files_in_folder(folder):
	files = []

	# r=root, d=directories, f = files
	for r, d, f in os.walk(folder):
		for file in f:
			files.append(file)

	return files


def list_files_in_folder(folder, extension):
	files = []

	# r=root, d=directories, f = files
	for r, d, f in os.walk(folder):
		for file in f:
			if file.endswith(extension):
				files.append( folder + "/" + file)

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


def merge_CSV(csv_file1, csv_file2, output_csv):
	all_runs = []

	# Open file 1
	reader1 = open_csv_file(csv_file1)
	row1 = next(reader1)
	all_runs.append(row1)

	# Open file 2
	reader2 = open_csv_file(csv_file2)
	row2 = next(reader2)
	all_runs.append(row2)


	# Write out file
	write_csv(output_csv, all_runs)


def append_csv(csv_file, row):

	# Open file
	reader = open_csv_file(csv_file)

	# Append row
	all_runs = []
	row = next(reader)
	all_runs.append(row)

	# Write out file
	write_csv(csv_file, all_runs)


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


def extract_csv_cell(filename, delimiter_char, col, row):
	if ( not os.path.isfile(filename) ):
		print( filename + " does not exist! ")
		return None


	with open(filename) as csv_file:
		try:
			data = [row for row in csv.reader(csv_file, delimiter=delimiter_char)]
			return data[col][row]

		except ValueError as e:
			print ("CSV file " + filename + " is invalid! " + str(e))
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
			print ("Json file " + filename + " is invalid! " + str(e))
			return None
