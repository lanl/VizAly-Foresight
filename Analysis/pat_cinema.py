#! /usr/bin/env python
""" This executable reads a JSON file and appends the metrics file with paths to images.
The images are plotted in this executable.
"""

import argparse
import csv
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy
import os
from pat import cinema
from pat import file_utilities as futils
from pat import plot_utilities as putils
from pat import Job as j

class PATCinema(cinema.CinemaWorkflow):

    def prepare_cinema(self):

        # open CSV file
        metrics_csv = self.json_data["project-home"] +  self.json_data["wflow-path"] + "/cbench/" + self.json_data["cbench"]["output"]["metrics-file"] + ".csv"

        # get list of all images
        image_data = {}
        for sec in self.json_data["pat"]["analysis"]:
            col_name = sec["output-column"]
            prefix = sec["output-prefix"]
            path = sec["path"]
            if col_name not in image_data.keys():
                image_data[col_name] = {}
            image_data[col_name][prefix] = path
        image_columns = list(image_data.keys())
        image_columns.sort()

        # loop over lines in metrics file and append image column and the image paths 
        all_runs = []
        with open(metrics_csv,"r") as csvinput:
            reader = csv.reader(csvinput)

            # modify Cinema files
            row = next(reader)
            for col_name in image_columns:
                row.append(col_name)
            all_runs.append(row)

            # loop over rows
            for row in reader:
                prefix = row[1].strip()

                # loop over image types
                for col_name in image_columns:
                    print(prefix)

                    # get data file
                    # if no data file recorded then continue
                    if prefix in image_data[col_name].keys():
                        data_file = image_data[col_name][prefix]
                    else:
                        row.append("")
                        continue

                    # see if data file exists
                    if not os.path.exists(data_file):
                        row.append("")
                        continue

                    # get original data
                    orig = numpy.loadtxt(data_file.replace(prefix, "orig"))

                    # plot
                    data = numpy.loadtxt(data_file)
                    plt.plot(data[:, 0], data[:, 1] / orig[:, 1])
                    output_file = col_name + "_" + prefix + ".png"
                    plt.savefig(output_file)
                    plt.close()

                    # append value to row
                    row.append(output_file)

                # store appended row for output
                all_runs.append(row)

            # write data CSV file
            futils.write_csv("data.csv", all_runs)

# parse Input
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--input-file", required=True)
parser.add_argument("--output-file", default="results.cdb")
opts = parser.parse_args()

# load data
cinema = PATCinema(opts.input_file)

# create directory
if not os.path.exists(opts.output_file):
     os.mkdir(opts.output_file)
os.chdir(opts.output_file)

# generate plots
cinema.prepare_cinema()
