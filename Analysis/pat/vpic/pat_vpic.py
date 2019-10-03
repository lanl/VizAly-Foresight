#! /usr/bin/env python
""" This script generates a workflow that runs CBench and VPIC data analysis executables.
"""

import argparse
import sys
import os

from pat.utils import file_utilities as futils
from pat.utils import job as j
from pat.utils import workflow

class VPICWorkflow(workflow.Workflow):

    def add_analysis_jobs(self):

        # get CBench job which is parent to all jobs in this function
        cbench_job = self.jobs[0]

        # get environment script
        environment = self.environment_from_json_data()

        # get spectrum information from configuration file
        spectrum_section = "spectrum"
        spectrum_config = self.configuration_from_json_data(spectrum_section)
        spectrum_exe = self.json_data["pat"]["analysis-tool"]["analytics"] \
                                     [spectrum_section]["path"]

        # create jobs to run analysis jobs on each output file from CBench
        for path in self.json_data["pat"]["input-files"]:
            print("Creating analysis jobs for", path)
            prefix = path["output-prefix"]
            cbench_path = path["path"]
            timestep = cbench_path.split(".")[-1]

            # create job for power spectrum
            spectrum_file = self.workflow_dir + "/" + spectrum_section + "/{}_spectrum.csv".format(prefix)
            spectrum_job = j.Job(name="{}_{}".format(prefix, spectrum_section),
                                 execute_dir=spectrum_section,
                                 executable=spectrum_exe,
                                 arguments=["--input-file", cbench_path, "--output-file", spectrum_file],
                                 configurations=spectrum_config,
                                 environment=environment)

            # make dependent on CBench job and add to workflow
            spectrum_job.add_parents(cbench_job)
            self.add_job(spectrum_job)

            # track output files
            self.json_data["pat"]["analysis"].append({"output-column" : "FILE_Spectrum",
                                                      "output-prefix" : prefix, "path" : spectrum_file})

    def add_cinema_plotting_jobs(self):

        # get environment for Cinema job
        if "evn_path" in self.json_data["cinema-plots"]:
                environment = self.json_data["foresight-home"] + self.json_data["cinema-plots"]["evn_path"]
        else:
                environment = None

        # get configuration for Cinema job
        if "configuration" in self.json_data["cinema-plots"]:
                configurations = list(sum(self.json_data["cinema-plots"]["configuration"].items(), ()))
        else:
                configurations = None

        # get executable for Cinema job
        cinema_exe = self.json_data["cinema-plots"]["path"]

        # create job for Cinema
        cinema_job = j.Job(name="cinema",
                           execute_dir="cinema",
                           executable=cinema_exe,
                           arguments=["--input-file", self.json_path],
                           configurations=configurations,
                           environment=environment)

        # make dependent on all previous workflow jobs and add to workflow
        cinema_job.add_parents(*self.jobs)
        self.add_job(cinema_job)

# parse command line
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--input-file", required=True)
parser.add_argument("--submit", action="store_true")
opts = parser.parse_args()

# read input JSON file
wflow_data = futils.read_json(opts.input_file)

# create Workflow instance
wflow_dir = wflow_data["project-home"]
wflow = VPICWorkflow("wflow", wflow_data, workflow_dir=wflow_dir)

# make directories
if not os.path.exists(wflow_dir + "/cbench"):
    os.makedirs(wflow_dir + "/cbench")
if not os.path.exists(wflow_dir + "/spectrum"):
    os.makedirs(wflow_dir + "/spectrum")
if not os.path.exists(wflow_dir + "/cinema"):
    os.makedirs(wflow_dir + "/cinema")

# add jobs to workflow
wflow.add_cbench_job()
wflow.add_analysis_jobs()
wflow.add_cinema_plotting_jobs()

# write submission script
wflow.write_submit()

# submit workflow
if opts.submit:
    wflow.submit()
