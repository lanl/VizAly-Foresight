#! /usr/bin/env python
""" Creates a workflow for Slurm clusters to compress and analyze cosmological data.
"""

import argparse
import configparser
import itertools
import json
import numpy
import os

class Workflow(object):
    """ Class that describes workflow for Slurm scheduler. Note that jobs ordering appended to workflow
    has meaning in this simple API.
    """

    def __init__(self, name, workflow_dir=""):

        # store meta-data about the workflow
        self.name = name

        # store a list of jobs in workflow
        self.jobs = []

        # location to store workflow files
        self.workflow_dir = workflow_dir

        # attributes for workflow construction
        self.submit_file = None

    def add_job(self, job, dependencies=None):
        """ Adds a job to the workflow.
        """
        self.jobs.append(job)

    def write(self):
         """ Writes Slurm workflow files.
         """
    
         # create workflow directory
         if len(self.workflow_dir):
             if not os.path.exists(self.workflow_dir):
                 os.mkdir(self.workflow_dir)
    
         # create submission script
         self.submit_path = self.name + ".sh"
         with open(self.submit_path, "w") as fp:
             fp.write("#! /bin/bash\n")
    
         # loop over each job
         for i, job in enumerate(self.jobs):
    
             # get a unique name and index
             job.name = job.name if job.name != None else "job_{}".format(i)
             job._idx = i
    
             # figure out dependencies job indices
             parent_idxs = [j._idx for j in job._parents]
             if parent_idxs == []:
                 depends_str = ""
             else:
                 depends_str = "--dependency=afterok:" + ":".join(["$jid{}".format(j) for j in parent_idxs])
    
             # create workflow directory
             if not os.path.exists(job.execute_dir):
                 os.mkdir(job.execute_dir)
    
             # write wrapper script for job
             path = self.workflow_dir + "/" + job.execute_dir + "/" + job.name + ".sh"
             with open(path, "w") as fp:
                 fp.write("#! /bin/bash\n")
                 for key, val in zip(job.configurations[::2], job.configurations[1::2]):
                     fp.write("#SBATCH --{}={}\n".format(key, val))
                 fp.write("mkdir -p {}/{}\n".format(self.workflow_dir, job.execute_dir))
                 fp.write("cd {}/{}\n".format(self.workflow_dir, job.execute_dir))
                 if job.environment != None:
                     fp.write("source {}\n".format(job.environment))
                 fp.write(job.executable + " " + " ".join(map(str, job.arguments)) + "\n")
    
             # append job to controller file
             slurm_out_path = self.workflow_dir + "/" + job.execute_dir + "/" + job.name + ".slurm.out"
             with open(self.submit_path, "a") as fp:
                 fp.write("\n# {}\n".format(job.name))
                 fp.write("jid{}=$(sbatch {} --output {} {})\n".format(job._idx, depends_str, slurm_out_path, path))
                 fp.write("jid{}=$(echo $jid{} | rev | cut -f 1 -d ' ' | rev)".format(job._idx, job._idx))
    
    def submit(self):
        """ Submits Slurm workflow.
        """
        os.system("bash {}".format(self.submit_path))

class Job(object):
    """ Class that describes a single job in a workflow.
    """

    def __init__(self, executable, execute_dir=".", arguments=None, configurations=[], name=None,
                 environment=None):

        # store meta-data about the job
        self.name = name
        self._idx = None
        self.configurations = configurations if configurations != None else []
        self.environment = environment

        # store command
        self.execute_dir = execute_dir
        self.executable = executable
        self.arguments = arguments if arguments != None else []

        # store dependencies
        self._parents = []
        self._childs = []

    def __repr__(self):
        return str(self.name)

    def add_parents(self, *jobs):
        self._parents += jobs
        for job in jobs:
            job._childs.append(self)

class Config(configparser.ConfigParser):
    """ Class for reading configuration file.
    """

    def __init__(self):
        super().__init__()

    def geteval(self, option, key):
        """ Does `eval` on getter.
        """
        return eval(self.get(option, key))

    def arguments_from_section(self, section):
        """ Returns list of arguments for section.
        """
        if cp.has_section(section):
            its = cp.items(section)
            its = [("--" + k, v) for k, v in its]
            return list(itertools.chain(its))
        return None

    def configuration_from_section(self, section):
        """ Returns list of configuration settings for section.
        """
        section += "-configuration"
        if cp.has_section(section):
            return list(itertools.chain(*cp.items(section)))
        return None

    def environment_from_section(self, section):
        """ Returns enviroment file for section.
        """
        if cp.has_section("environments"):
            if cp.has_option("enviroments", section):
                return cp.get("environments", section)
        return None

# parse command line
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--name", default="workflow_c", help="Name of the workflow.")
parser.add_argument("--config-file", default="workflow_c.ini", help="Path to configuration file.")
parser.add_argument("--submit", action="store_true", help="Submit workflow to Slurm.")
opts = parser.parse_args()

# read configuration file
cp = Config()
cp.readfp(open(opts.config_file))

# change directory
if not os.path.exists(opts.name):
    os.makedirs(opts.name)
else:
    raise OSError("The directory {} already exists! Aborting!".format(opts.name))
os.chdir(opts.name)
run_dir = os.getcwd()

# create output directories
cbench_dir = os.path.join(run_dir, "cbench")
halo_dir = os.path.join(run_dir, "halo")
spectra_dir = os.path.join(run_dir, "spectra")
output_dir = os.path.join(run_dir, "output")
os.makedirs(cbench_dir)
os.makedirs(halo_dir)
os.makedirs(spectra_dir)
os.makedirs(output_dir)

# create a workflow
wflow = Workflow(name=opts.name)

# template for building CBench JSON data
section = "cbench"
cbench_json_data = {
    "input" : {
        "filetype-comment" : "Type of file to load; HACC or NYX",
        "filetype" : cp.get(section, "file-type"),
        "filename-comment" : "Name of input file",
        "filename" : cp.get(section, "input-file"),
        "scalars-comment" : "Scalars to test",
        "scalars" : cp.geteval(section, "scalars"),
    },
    "output" : {
        "output-decompressed" : cp.getboolean(section, "output-decompressed"),
        "logfname-comment" : "Name of output log file",
        "metricsfname-comment" : "Name of file with output",
    },
    "compressor-comment" : "Compressors and parameters to test",
    "compressors" : [],
    "metrics-comment": "Metrics to report",
    "metrics": [],
}

# add list of metrics to compute with CBench JSON data
for metric in cp.geteval(section, "metrics"):
    entry = {"name" : metric}
    if cp.getboolean("cbench", "histogram") and metric == "absolute_error":
        entry["histogram"] = cp.geteval(section, "scalars")
    cbench_json_data["metrics"].append(entry)

# store files
cbench_files = []
halo_finder_files = []
spectra_files = []
metrics_files = []
plot_halo_dist_files = []
plot_spectra_files = []

# store compressor data
compressor_data = []
compressor_inputs = []

# loop over compressors
# first "compressor" is the input file
compressors = cp.items("compressors")
compressors = [("original", None,)] + compressors
for i, (c_tag, c_name) in enumerate(compressors):

    # create a CBench job if doing compression
    if c_tag != "original":

        # get cartesian product of compressor settings to sweep
        keys = []
        vals = []
        for key, val in cp.items(c_tag):
            keys.append(key)
            val = eval(val)
            if isinstance(val, int) or isinstance(val, float):
                val = [val]
            vals.append(val)
        settings = list(itertools.product(*vals))
    
        # add compressors settings to JSON data
        cbench_json_data["compressors"] = []
        compressor_data.append({"name" : c_name})
        for i, setting in enumerate(settings):
            v_tag = ""
            for key, val in zip(keys, setting):
                v_tag += "{}:{}__".format(key, val)
                compressor_data[-1][key] = val
                if key not in compressor_inputs:
                    compressor_inputs.append(key)
            entry = {
                "name" : c_name,
                "output-prefix" : "out__{}__{}".format(c_tag, v_tag),
            }
            for j, val in enumerate(setting):
                entry[keys[j]] = val
            cbench_json_data["compressors"].append(entry)
    
        # set log file prefixes for CBench in JSON data
        section = "cbench"
        cbench_json_data["output"]["logfname"] = cp.get(section, "log-file") + "_{}".format(c_tag)
        cbench_json_data["output"]["metricsfname"] = cp.get(section, "metrics-file") + "_{}".format(c_tag)
        metrics_files.append(os.path.join(cbench_dir, cbench_json_data["output"]["metricsfname"] + ".csv"))

        # write CBENCH JSON data
        json_file = os.path.join(cbench_dir, "cbench_{}.json".format(c_tag))
        with open(json_file, "w") as fp:
            json.dump(cbench_json_data, fp, indent=4, sort_keys=True)
    
        # add a single CBench job to workflow for entire sweep
        cbench_job = Job(name="cbench_{}".format(c_tag),
                         execute_dir=cbench_dir,
                         executable=cp.get("executables", "mpirun"),
                         arguments=[cp.get("executables", section), json_file],
                         configurations=cp.configuration_from_section(section),
                         environment=cp.environment_from_section(section))
        wflow.add_job(cbench_job)

    # fill in JSON data if input file
    else:
        cbench_json_data["compressors"] = [{"name" : "original", "output-prefix" : "out__original__"}]

    # loop over each compressed file from CBench
    for i, _ in enumerate(cbench_json_data["compressors"]):

        # compressed file not explicitly set so construct it here
        # cut off timestep from path for halo finder executable
        if c_tag == "original":
            cbench_file = cp.get("cbench", "input-file")
            prefix = ".".join(cbench_file.split(".")[:-1])
        else:
            cbench_file = cbench_json_data["compressors"][i]["output-prefix"] + "__" + os.path.basename(cbench_json_data["input"]["filename"])
            prefix = cbench_dir + "/" + ".".join(cbench_file.split(".")[:-1])
        timestep = cbench_file.split(".")[:-1]
        cbench_files.append(cbench_file)

        # set paths for halo finder configuration and parameters files
        config_file = os.path.join(halo_dir, "halo_finder_{}_{}_config.txt".format(c_tag, i))
        parameters_file = os.path.join(halo_dir, "halo_finder_{}_{}_parameters.txt".format(c_tag, i))

        # write halo finder parameters file
        # specify location of parsed configuration file inside
        section = "halo-finder"
        os.system("sed \"s/^COSMOTOOLS_CONFIG.*/COSMOTOOLS_CONFIG .\/{}/\" {} > {}".format(os.path.basename(config_file),
                                                                                           cp.get(section, "parameters-file"),
                                                                                           parameters_file))

        # write halo finder configuration file
        # specify output prefix inside
        os.system("sed \"s/^BASE_OUTPUT_FILE_NAME.*/BASE_OUTPUT_FILE_NAME .\/{}/\" {} > {}".format(cbench_json_data["compressors"][i]["output-prefix"],
                                                                                                   cp.get(section, "config-file"),
                                                                                                   "tmp.out"))
        os.system("sed \"s/^ACCUMULATE_CORE_NAME.*/ACCUMULATE_CORE_NAME .\/{}/\" {} > {}".format(cbench_json_data["compressors"][i]["output-prefix"],
                                                                                                "tmp.out", config_file))

        # halo finder file not explicitly set so construct it here
        halo_finder_file = os.path.join(cbench_json_data["compressors"][i]["output-prefix"], "{}.fofproperties".format(timestep))
        halo_finder_files.append(halo_finder_file)

        # add halo finder job to workflow for compressed file
        # make dependent on CBench job
        halo_finder_job = Job(name="halo_finder_{}_{}".format(c_tag, i),
                              execute_dir=halo_dir,
                              executable=cp.get("executables", "mpirun"),
                              arguments=[cp.get("executables", section),
                                         "--config", config_file,
                                         "--timesteps", cp.get(section, "timesteps-file"),
                                         "--prefix", prefix,
                                         parameters_file],
                              configurations=cp.configuration_from_section(section),
                              environment=cp.environment_from_section(section))
        if c_tag != "original":
            halo_finder_job.add_parents(cbench_job)
        wflow.add_job(halo_finder_job)

        # spectra file not explicitly set so construct it here
        spectra_file = os.path.join(spectra_dir, "spectra_{}_{}.pk".format(c_tag, i))
        spectra_files.append(spectra_file)

        # add power spectra job to workflow for compressed file
        # make dependent on CBench job
        section = "power-spectrum"
        spectra_job = Job(name="spectra_{}_{}".format(c_tag, i),
                          execute_dir=spectra_dir,
                          executable=cp.get("executables", "mpirun"),
                          arguments=[cp.get("executables", section), cp.get(section, "parameters-file"),
                                     "-n", os.path.join(cbench_dir, cbench_file), spectra_file.rstrip(".pk")],
                                     configurations=cp.configuration_from_section(section),
                          environment=cp.environment_from_section(section))
        if c_tag != "original":
            spectra_job.add_parents(cbench_job)
        wflow.add_job(spectra_job)

        # add plot power spectra job to workflow for compressed file
        # make dependent on power spectra job
        section = "plot-halo-distribution"
        plot_halo_dist_files.append(output_dir + "/plot_halo_dist_{}_{}.png".format(c_tag, i))
        plot_halo_dist_job = Job(name="plot_halo_dist_{}_{}".format(c_tag, i),
                                 execute_dir=output_dir,
                                 executable=cp.get("executables", section),
                                 arguments=["--input-files", halo_finder_files[-1], "--reference-file", halo_finder_files[0],
                                            "--output-file", plot_halo_dist_files[-1]] + cp.arguments_from_section(section),
                                 configurations=cp.configuration_from_section(section),
                                 environment=cp.environment_from_section(section))
        plot_halo_dist_job.add_parents(halo_finder_job)
        wflow.add_job(plot_halo_dist_job)

        # add plot power spectra job to workflow for compressed file
        # make dependent on power spectra job
        section = "plot-power-spectrum"
        plot_spectra_files.append(output_dir + "/plot_spectra_{}_{}.png".format(c_tag, i))
        plot_spectra_job = Job(name="plot_spectra_{}_{}".format(c_tag, i),
                               execute_dir=output_dir,
                               executable=cp.get("executables", section),
                               arguments=["--input-files", spectra_files[-1], "--reference-file", spectra_files[0],
                                          "--output-file", plot_spectra_files[-1]] + cp.arguments_from_section(section),
                               configurations=cp.configuration_from_section(section),
                               environment=cp.environment_from_section(section))
        plot_spectra_job.add_parents(spectra_job)
        wflow.add_job(plot_spectra_job)

# write workflow
wflow.write()

# write output manifest
with open(run_dir + "/{}.csv".format(opts.name), "w") as fp:
    fp.write(",".join(["compressor_name"] + compressor_inputs + ["metric_file", "halo_dist_file", "spectra_file"]) + "\n")
    lines = zip(metrics_files, plot_halo_dist_files, plot_spectra_files)
    for i, line in enumerate(lines):
        inputs = compressor_data[i]["name"] + ","
        for key in compressor_inputs:
            if key in compressor_data[i].keys():
                inputs += str(compressor_data[i][key]) + ","
            else:
                inputs += ","
        line = inputs + ",".join(line) + "\n"
        fp.write(line)

# submit
if opts.submit:
    wflow.submit()
