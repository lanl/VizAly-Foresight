#! /usr/bin/env python
""" This script generates a workflow that runs CBench and HACC data analysis executables.
"""

import argparse
import os
from pat import file_utilities as futils
from pat import Job as j
from pat import workflow

class HACCWorkflow(workflow.Workflow):

    def add_analysis_jobs(self):

        # get CBench job which is parent to all jobs in this function
        cbench_job = self.jobs[0]

        # get halo finder information from configuration file
        halo_section = "halo"
        halo_exe = self.json_data["simulation-analysis"]["analysis-tool"]["analytics"] \
                                 [halo_section]["path"]
        timesteps_path = self.json_data["simulation-analysis"]["analysis-tool"]["analytics"] \
                                       [halo_section]["timesteps-file"]
        config_path = self.json_data["simulation-analysis"]["analysis-tool"]["analytics"] \
                                    [halo_section]["config-file"]
        parameters_path = self.json_data["simulation-analysis"]["analysis-tool"]["analytics"] \
                                        [halo_section]["parameters-file"]

        # get halo finder run specifications from configuration file
        halo_config = self.configuration_from_json_data(halo_section)
        environment = self.environment_from_json_data()

        # get spectrum information from configuration file
        spectrum_section = "spectrum"
        spectrum_exe = self.json_data["simulation-analysis"]["analysis-tool"]["analytics"] \
                                     [spectrum_section]["path"]

        # get halo finder run specifications from configuration file
        spectrum_config = self.configuration_from_json_data(spectrum_section)

        # create job to run halo finder on each output file from CBench
        # create job to run 
        for path in self.json_data["simulation-analysis"]["input-files"]:
            print("Creating analysis jobs for", path)

            # create job for halo finder
            prefix = path["output-prefix"]
            halo_job = j.Job(name="{}_{}".format(prefix, halo_section),
                           execute_dir=halo_section,
                           executable=halo_exe,
                           arguments=["--config", config_path,
                                      "--timesteps", timesteps_path,
                                      "--prefix", prefix,
                                      parameters_path],
                           configurations=halo_config,
                           environment=environment)

            # make dependent on CBench job and add to workflow
            halo_job.add_parents(cbench_job)
            self.add_job(halo_job)

            # create job for power spectrum
            spectrum_job = j.Job(name="{}_{}".format(prefix, spectrum_section),
                                 execute_dir=spectrum_section,
                                 executable=spectrum_exe,
                                 arguments=["-n", "FILE", 499],
                                 configurations=spectrum_config,
                                 environment=environment)

            # make dependent on CBench job and add to workflow
            spectrum_job.add_parents(cbench_job)
            self.add_job(spectrum_job)

    def add_plotting_jobs(self):
        print("PLOTTING JOBS")

        halo_images = []
        spectrum_images = []

        #create_CinemaDB()

# parse command li9ne
parser = argparse.ArgumentParser()
parser.add_argument("--input-file")
parser.add_argument("--submit", action="store_true")
opts = parser.parse_args()

# read input JSON file
wflow_data = futils.read_json(opts.input_file)

# create Workflow instance
wflow_dir = wflow_data["project-home"]
wflow = HACCWorkflow("wflow", wflow_data, workflow_dir=wflow_dir)

# add jobs to workflow
wflow.add_cbench_job()
wflow.add_analysis_jobs()
wflow.add_plotting_jobs()

# write submission script
wflow.write_submit()

# submit workflow
if opts.submit:
    wflow.submit()
