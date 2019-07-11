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
        halo_config = self.configuration_from_json_data(halo_section)
        halo_exe = self.json_data["pat"]["analysis-tool"]["analytics"] \
                                 [halo_section]["path"]
        timesteps_path = self.json_data["pat"]["analysis-tool"]["analytics"] \
                                       [halo_section]["timesteps-file"]
        config_path = self.json_data["pat"]["analysis-tool"]["analytics"] \
                                    [halo_section]["config-file"]
        parameters_path = self.json_data["pat"]["analysis-tool"]["analytics"] \
                                        [halo_section]["parameters-file"]

        # get halo distribution information from configuraiton file
        halo_query_section = "halo_query"
        halo_query_config = self.configuration_from_json_data(halo_query_section)
        halo_query_exe = self.json_data["pat"]["analysis-tool"]["analytics"] \
                                 [halo_query_section]["path"]

        # get spectrum information from configuration file
        spectrum_section = "spectrum"
        spectrum_config = self.configuration_from_json_data(spectrum_section)
        spectrum_exe = self.json_data["pat"]["analysis-tool"]["analytics"] \
                                     [spectrum_section]["path"]
        spectrum_config_path = self.json_data["pat"]["analysis-tool"]["analytics"] \
                                             [spectrum_section]["config-file"]

        # create jobs to run analysis jobs on each output file from CBench
        for path in self.json_data["pat"]["input-files"]:
            print("Creating analysis jobs for", path)
            prefix = path["output-prefix"]
            cbench_path = path["path"]
            timestep = cbench_path.split(".")[-1]

            # write halo finder parameters file
            # specify location of parsed configuration file inside
            new_parameters_path = halo_section + "/" + prefix + "_halo_params.txt"
            new_config_path = halo_section + "/" + prefix + "_halo_config.txt"
            os.system("sed \"s/^COSMOTOOLS_CONFIG.*/COSMOTOOLS_CONFIG .\/{}/\" {} > {}".format(
                                    os.path.basename(new_config_path),
                                    parameters_path,
                                    self.workflow_dir + "/" + new_parameters_path))

            # write halo finder configuration file
            # specify output prefix inside
            tmp_path = "tmp.out"
            os.system("sed \"s/^BASE_OUTPUT_FILE_NAME.*/BASE_OUTPUT_FILE_NAME .\/{}/\" {} > {}".format(
                                    prefix, config_path, tmp_path))
            os.system("sed \"s/^ACCUMULATE_CORE_NAME.*/ACCUMULATE_CORE_NAME .\/{}/\" {} > {}".format(
                                    prefix, tmp_path, self.workflow_dir + "/" + new_config_path))
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

            # predict halo finder file
            halo_file = self.workflow_dir + "/" + halo_section + "/" + prefix + "-" + timestep + ".fofproperties"

            # loop over different halo SQL queries
            # create job for getting halo distribution
            sec = self.json_data["pat"]["analysis-tool"]["analytics"][halo_query_section]
            for i, _ in enumerate(sec["query"]):
                halo_query_file = self.workflow_dir + "/" + halo_query_section + "/halo_query_{}_{}.csv".format(prefix, i)
                args = ["--input-file", halo_file,
                        "--output-file", halo_query_file,
                        "--query", "\"{}\"".format(sec["query"][i]),
                        "--xlabel", "\"{}\"".format(sec["xlabel"][i]),
                        "--ylabel", "\"{}\"".format(sec["ylabel"][i])]
                if "xlim" in sec.keys():
                    args += ["--xlim", " ".join(map(str, sec["xlim"][i]))]
                if "log-bins" in sec.keys():
                    if sec["log-bins"][i]:
                        args += ["--log-bins"]
                halo_query_job = j.Job(name="{}_{}_{}".format(prefix, halo_query_section, i),
                                       execute_dir=halo_query_section,
                                       executable=halo_query_exe,
                                       arguments=args,
                                       configurations=halo_query_config,
                                       environment=environment)
                self.json_data["pat"]["analysis"].append({"output-column" : "FILE_Halo_Distribution_{}".format(i),
                                                          "output-prefix" : prefix, "path" : halo_query_file})   
 
                # make dependent on halo finder job
                halo_query_job.add_parents(halo_job)
                self.add_job(halo_query_job)

            # create job for power spectrum
            spectrum_job = j.Job(name="{}_{}".format(prefix, spectrum_section),
                                 execute_dir=spectrum_section,
                                 executable=spectrum_exe,
                                 arguments=[spectrum_config_path, "-n", cbench_path, prefix + "_spectrum", timestep],
                                 configurations=spectrum_config,
                                 environment=environment)

            # make dependent on CBench job and add to workflow
            spectrum_job.add_parents(cbench_job)
            self.add_job(spectrum_job)

            # predict and track output files
            for ext in ["", ".rsd.0", ".rsd.1", ".rsd.2"]:
                spectrum_file = self.workflow_dir + "/" + spectrum_section + "/{}_spectrum.pk{}".format(prefix, ext)
                self.json_data["pat"]["analysis"].append({"output-column" : "FILE_Spectrum{}".format(ext),
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
                           arguments=["--input-file", self.json_file],
                           configurations=configurations,
                           environment=environment)

        # make dependent on all previous workflow jobs and add to workflow
        cinema_job.add_parents(*self.jobs)
        self.add_job(cinema_job)

# parse command line
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
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
os.makedirs(wflow_dir + "/cinema")

# add jobs to workflow
wflow.add_cbench_job()
wflow.add_analysis_jobs()

# write submission script
wflow.write_submit()

# submit workflow
if opts.submit:
    wflow.submit()
