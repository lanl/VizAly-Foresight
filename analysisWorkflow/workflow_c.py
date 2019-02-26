#! /usr/bin/env python

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
        """ Does eval on getter.
        """
        return eval(self.get(option, key))

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument("--name", default="workflow_c")
parser.add_argument("--config-file", default="workflow_c.ini")
parser.add_argument("--submit", action="store_true")
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
os.makedirs(cbench_dir)
os.makedirs(halo_dir)
os.makedirs(spectra_dir)

# create a workflow
wflow = Workflow(name=opts.name)

# template for building CBench JSON file
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
for metric in cp.geteval(section, "metrics"):
    entry = {"name" : metric}
    if cp.getboolean("cbench", "histogram") and metric == "absolute_error":
        entry["histogram"] = cp.geteval(section, "scalars")
    cbench_json_data["metrics"].append(entry)

# loop over compressor settings
for c_tag, c_name in cp.items("compressors"):

    # get cartesian product of compressors settings
    # specify output prefix
    cbench_json_data["compressors"] = []
    keys = []
    vals = []
    for key, val in cp.items(c_tag):
        keys.append(key)
        val = eval(val)
        if isinstance(val, int) or isinstance(val, float):
            val = [val]
        vals.append(val)
    settings = list(itertools.product(*vals))
    for i, setting in enumerate(settings):
        entry = {
            "name" : c_name,
            "output-prefix" : "__{}__{}".format(c_tag, i),
        }
        for i, val in enumerate(setting):
            entry[keys[i]] = val
        cbench_json_data["compressors"].append(entry)

    # set log files for CBench
    section = "cbench"
    cbench_json_data["output"]["logfname"] = cp.get(section, "log-file") + "_{}".format(c_tag)
    cbench_json_data["output"]["metricsfname"] = cp.get(section, "metrics-file") + "_{}".format(c_tag)

    # write JSON data
    json_file = os.path.join(cbench_dir, "cbench_{}.json".format(c_tag))
    with open(json_file, "w") as fp:
        json.dump(cbench_json_data, fp, indent=4, sort_keys=True)

    # add CBench job to workflow
    cbench_job = Job(name="cbench_{}".format(c_tag),
                     execute_dir=cbench_dir,
                     executable=cp.get("executables", "mpirun"),
                     arguments=[cp.get("executables", section), json_file],
                     configurations=list(itertools.chain(*cp.items("{}-configuration".format(section)))),
                     environment=cp.get(section, "environment-file") if cp.get(section, "environment-file") else None)
    wflow.add_job(cbench_job)

    # loop over each compressed file from CBench
    for i, setting in enumerate(settings):

        # get CBench output path
        cbench_file = cbench_json_data["compressors"][i]["output-prefix"] + "__" + os.path.basename(cbench_json_data["input"]["filename"])

        # cut off timestep from CBench output path for halo finder executable
        prefix = cbench_dir + "/" + ".".join(cbench_file.split(".")[:-1])

        # add halo finder job to workflow
        # make dependent on CBench job
        section = "halo-finder"
        args = [cp.get("executables", section),
                "--config", cp.get(section, "config-file"),
                "--timesteps", cp.get(section, "timesteps-file"),
                "--prefix", prefix,
                cp.get(section, "parameters-file")]
        halo_finder_job = Job(name="halo_finder_{}_{}".format(c_tag, i),
                              execute_dir=halo_dir,
                              executable=cp.get("executables", "mpirun"),
                              arguments=args,
                              configurations=list(itertools.chain(*cp.items("{}-configuration".format(section)))),
                              environment=cp.get(section, "environment-file") if cp.get(section, "environment-file") else None)
        halo_finder_job.add_parents(cbench_job)
        wflow.add_job(halo_finder_job)

        # add power spectra job to workflow
        section = "power-spectrum"        
        spectra_job = Job(name="spectra_{}_{}".format(c_tag, i),
                          execute_dir=spectra_dir,
                          executable=cp.get("executables", "mpirun"),
                          arguments=[cp.get("executables", section), cp.get(section, "parameters-file"),
                                     "-n", os.path.join(cbench_dir, cbench_file),
                                     os.path.join(spectra_dir, "spectra_{}_{}".format(c_tag, i))],
                              configurations=list(itertools.chain(*cp.items("{}-configuration".format(section)))),
                              environment=cp.get(section, "environment-file") if cp.get(section, "environment-file") else None)
        spectra_job.add_parents(cbench_job)
        wflow.add_job(spectra_job)

# write workflow
wflow.write()

# submit
if opts.submit:
    wflow.submit()
