#!/usr/bin/python

import sys, os, json, csv

import file_utilities as f
import job as j


class Analysis_Workflow:
	sim = ""

	def __init__(self, sim_name):
		self.sim = sim_name

	def init_analysis(self, num_input_params, json_file):
		f.validate_num_args(num_input_params, 2)
		self.json_data =  f.read_json(json_file)

	def set_analysis(self, json_file):
		self.json_data =  f.read_json(json_file)

	def launch_sim_analysis(self, func_name):
		return func_name(self.json_data)
