#!/usr/bin/python

import sys
import os
import json
import csv
from collections import OrderedDict

import Workflow.draw.utils as utils

class CinemaCreator:

	def __init__(self, json_file):
		self.json_file = json_file
		print(self.json_file)
		self.json_data = utils.read_json(self.json_file)


	def create_cinema_database(self, cinema_database, csv_file, image_files):
		#create cdb file
		utils.create_folder(cinema_database)

		# Copy img files to cinema
		for img in image_files:
			cmd = "cp " + img + " " + cinema_database
			os.system(cmd)

		# Copy data to cinema
		cmd = "cp " + csv_file + " " + cinema_database + "/data.csv"
		os.system(cmd)

		# Copy JSON to cinema for provenance
		cmd = "cp " + self.json_file + " " + cinema_database + "/wflow.json"
		os.system(cmd)

		print("Create cinema database " + cinema_database)


	def create_cinema(self):
		self.prepare_cinema()

		# Create cinema database
		output_path = self.json_data['project-home'] + self.json_data['wflow-path'] + "/"
		cinema_database = output_path + self.json_data['visualize']['cinema']['name'] + ".cdb"
		img_files = utils.list_files_in_folder( output_path + "plots/" , ".png")
		cinema_csv = output_path + "cinema/" + "data.csv"
		print(output_path)
		print(img_files)

		self.create_cinema_database(cinema_database, cinema_csv, img_files)


	def prepare_cinema(self):
		# Adds plotting jobs to workflow that produce final products.
		raise NotImplementedError("Implement the `prepare_cinema` function to your workflow!")