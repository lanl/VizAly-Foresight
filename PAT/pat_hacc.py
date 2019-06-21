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

        # get environment script
        environment = self.environment_from_json_data()

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

        # get spectrum information from configuration file
        spectrum_section = "spectrum"
        spectrum_exe = self.json_data["simulation-analysis"]["analysis-tool"]["analytics"] \
                                     [spectrum_section]["path"]
        spectrum_config_path = self.json_data["simulation-analysis"]["analysis-tool"]["analytics"] \
                                             [spectrum_section]["config-file"]

        # get halo finder run specifications from configuration file
        spectrum_config = self.configuration_from_json_data(spectrum_section)

        # create job to run halo finder on each output file from CBench
        # create job to run 
        for path in self.json_data["simulation-analysis"]["input-files"]:
            print("Creating analysis jobs for", path)
            prefix = path["output-prefix"]
            cbench_path = path["path"]
            timestep = 499

            # write halo finder parameters file
            # specify location of parsed configuration file inside
            new_parameters_path = halo_section + "/" + prefix + "_halo_params.txt"
            new_config_path = halo_section + "/" + prefix + "_halo_config.txt"
            os.system("sed \"s/^COSMOTOOLS_CONFIG.*/COSMOTOOLS_CONFIG .\/{}/\" {} > {}".format(os.path.basename(new_config_path),
                                                                                               parameters_path,
                                                                                               self.workflow_dir + "/" + new_parameters_path))

            # write halo finder configuration file
            # specify output prefix inside
            tmp_path = "tmp.out"
            os.system("sed \"s/^BASE_OUTPUT_FILE_NAME.*/BASE_OUTPUT_FILE_NAME .\/{}/\" {} > {}".format(prefix,
                                                                                                       config_path,
                                                                                                       tmp_path))
            os.system("sed \"s/^ACCUMULATE_CORE_NAME.*/ACCUMULATE_CORE_NAME .\/{}/\" {} > {}".format(prefix,
                                                                                                     tmp_path,
                                                                                                     self.workflow_dir + "/" + new_config_path))
            os.remove(tmp_path)

            # create job for halo finder
            halo_job = j.Job(name="{}_{}".format(prefix, halo_section),
                           execute_dir=halo_section,
                           executable=halo_exe,
                           arguments=["--config", os.path.basename(new_config_path),
                                      "--timesteps", timesteps_path,
                                      "--prefix", cbench_path[:-len(str(timestep)) - 1],
                                      os.path.basename(new_parameters_path)],
                           configurations=halo_config,
                           environment=environment)

            # make dependent on CBench job and add to workflow
            halo_job.add_parents(cbench_job)
            self.add_job(halo_job)

            # create job for power spectrum
            spectrum_job = j.Job(name="{}_{}".format(prefix, spectrum_section),
                                 execute_dir=spectrum_section,
                                 executable=spectrum_exe,
                                 arguments=[spectrum_config_path, "-n", cbench_path, prefix + "spectrum", timestep],
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

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument("--input-file")
parser.add_argument("--submit", action="store_true")
opts = parser.parse_args()

# read input JSON file
wflow_data = futils.read_json(opts.input_file)

# create Workflow instance
wflow_dir = wflow_data["project-home"]
wflow = HACCWorkflow("wflow", wflow_data, workflow_dir=wflow_dir)

# make directories
os.makedirs(wflow_dir + "/cbench")
os.makedirs(wflow_dir + "/halo")
os.makedirs(wflow_dir + "/spectrum")

# add jobs to workflow
wflow.add_cbench_job()
wflow.add_analysis_jobs()
wflow.add_plotting_jobs()

# write submission script
wflow.write_submit()

# submit workflow
if opts.submit:
    wflow.submit()
