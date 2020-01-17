#!/usr/bin/python
import sys
import os
import json
import csv

class Job(object):
	""" Class that describes a single job in a workflow.
	"""

	def __init__(self, executable, execute_dir=".", arguments=None, configurations=[], name=None, environment=None):

		# store meta-data about the job
		self.name = name
		self._idx = None
		self.configurations = configurations if configurations != None else []
		self.environment = environment

		# store command
		self.execute_dir = execute_dir
		self.executable = executable
		print("args", arguments)
		self.arguments = arguments if arguments != None else []
		self.commands = []

		# store wflow.json
		self._parents = []
		self._childs = []


	def __repr__(self):
		return str(self.name)


	def add_command(self, command):
		self.commands.append(command)


	def add_parents(self, *jobs):
		self._parents += jobs
		for job in jobs:
			job._childs.append(self)