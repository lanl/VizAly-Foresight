#!/usr/bin/python

import sys, os, json, csv
from collections import OrderedDict

import file_utilities as futils
import plot_utilities as putils



class CinemaWorkflow():

	def __init__(self, json_file):
		self.json_file = json_file
		print self.json_file
		self.json_data = futils.read_json(self.json_file)


	def create_cinema_database(self, cinema_database, csv_file, image_files):
		#create cdb file
		futils.create_folder(cinema_database)

		# Copy files to cinema
		for img in image_files:
			cmd = "cp " + img + " " + cinema_database
			os.system(cmd)

		# Copy data to cinema
		cmd = "cp " + csv_file + " " + cinema_database + "/data.csv"
		os.system(cmd)

		print "Create cinema database " + cinema_database


	def create_plots(self):
		pass


	def prepare_cinema(self):
		pass


	def create_cinema(self):
		self.prepare_cinema()

		# Create cinema database
		cinema_database = self.json_data['simulation-analysis']['analysis-folder'] + ".cdb"
		img_files = futils.list_files_in_folder( self.json_data['project-home'] + "/cinema/" , ".png")
		cinema_csv = self.json_data["project-home"] + "/cbench/" + "data.csv"

		self.create_cinema_database(cinema_database, cinema_csv, img_files)


	



