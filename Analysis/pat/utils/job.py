#!/usr/bin/python
import sys
import os
import json
import csv


class Job(object):
    """ Class that describes a single job in a workflow.
    """

    def __init__(self, executable, name=None, execute_dir=".", arguments=None, configurations=[], 
                 environment=None, job_type=None, commands=[]):

        # store meta-data about the job
        self.name = name
        self.job_type = job_type
        self._idx = None
        self.configurations = configurations if configurations != None else []
        self.environment = environment

        # store command
        self.execute_dir = execute_dir
        self.executable = executable
        self.arguments = arguments if arguments != None else []

        # specific instructions like create a folder
        self.commands = []
        self.post_commands = []

        # store wflow.json
        self._parents = []
        self._childs = []


    def __repr__(self):
        return str(self.name)


    def add_command(self, command):
        self.commands.append(command)


    def add_post_command(self, command):
        self.post_commands.append(command)


    def add_parent(self, job):
        self._parents.append(job)
        job._childs.append(self)
        

    def add_parents(self, jobs):
        for job in jobs:
            self._parents.append(job)

        for job in jobs:
            job._childs.append(self)
