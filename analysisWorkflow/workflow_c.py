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

    def __init__(self, executable, execute_dir=".", arguments=None, configurations=[], name=None):

        # store meta-data about the job
        self.name = name
        self._idx = None
        self.configurations = configurations if configurations != None else []

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

# template for building CBench JSON file
cbench_json_data = {
    "input" : {
        "filetype-comment" : "Type of file to load; HACC or NYX",
        "filetype" : None,
        "filename-comment" : "Name of input file",
        "filename" : None,
        "scalars-comment" : "Scalars to test",
        "scalars" : None,
    },
    "output" : {
        "output-decompressed" : None,
        "logfname-comment" : "Name of output log file",
        "logfname" : None,
        "metricsfname-comment" : "Name of file with output",
        "metricsfname" : None,
    },
    "compressor-comment" : "Compressors and parameters to test",
    "compressors" : [],
    "metrics-comment": "Metrics to report",
    "metrics": [],
}

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument("--name", default="workflow_c")
parser.add_argument("--config-file", default="workflow_c.ini")
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
os.makedirs(cbench_dir)

# create a workflow
wflow = Workflow(name=opts.name)

# fill in CBench JSON data except for list of compressors and settings
section = "cbench"
cbench_json_data["input"]["filetype"] = cp.get(section, "file-type")
cbench_json_data["input"]["filename"] = cp.get(section, "input-file")
cbench_json_data["input"]["scalars"] = cp.geteval(section, "scalars")
cbench_json_data["output"]["output-decompressed"] = cp.getboolean(section, "output-decompressed")
cbench_json_data["output"]["logfname"] = cp.get(section, "log-file")
cbench_json_data["output"]["metricsfname"] = cp.get(section, "metrics-file")
for metric in cp.geteval(section, "metrics"):
    entry = {"name" : metric}
    if cp.getboolean("cbench", "histogram") and metric == "absolute_error":
        entry["histogram"] = cp.geteval(section, "scalars")
    cbench_json_data["metrics"].append(entry)

# get CBench settings

# loop over compressor settings
for c_tag, c_name in cp.items("compressors"):

    # clear compressors
    cbench_json_data["compressors"] = []

    # get cartesian product of compressors settings
    keys = []
    vals = []
    for key, val in cp.items(c_tag):
        keys.append(key)
        val = eval(val)
        if isinstance(val, int) or isinstance(val, float):
            val = [val]
        vals.append(val)
    settings = list(itertools.product(*vals))
    for setting in settings:
        entry = {"name" : c_name}
        for i, val in enumerate(setting):
            entry[keys[i]] = val
        cbench_json_data["compressors"].append(entry)

    # write JSON data
    json_file = os.path.join(cbench_dir, "cbench_{}.json".format(c_tag))
    with open(json_file, "w") as fp:
        json.dump(cbench_json_data, fp, indent=4, sort_keys=True)

    # add CBench job to workflow
    cbench_job = Job(name="cbench_{}".format(c_tag),
                     execute_dir=cbench_dir,
                     executable=cp.get("executables", "cbench"),
                     arguments=[json_file],
                     configurations=list(itertools.chain(*cp.items("cbench-configuration"))))
    wflow.add_job(cbench_job)


# run halo finder

# run power spectra

# concatenate in Cinema

# write workflow
wflow.write()
