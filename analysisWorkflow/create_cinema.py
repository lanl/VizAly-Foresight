#! /usr/bin/env python
""" Creates a Cinema database from workflow outputs. Plots are created in serial.
"""

import argparse
import os
import pandas
import subprocess

def _external_call(cmd, debug=False):
    """ This function makes external calls.
    """
    cmd = list(map(str, cmd))
    if debug:
        print(" ".join(cmd))
    else:
        cmd += [">", "/dev/null", "2>&1"]
    os.system(" ".join(cmd))

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument("--input-file")
parser.add_argument("--output-file", default="out.cdb")
parser.add_argument("--reference-row", default="Original")
parser.add_argument("--debug", action="store_true")
opts = parser.parse_args()

# read manifest
data = pandas.read_csv(opts.input_file)

# make output file
if not os.path.exists(opts.output_file):
    os.makedirs(opts.output_file)

# get CSV columns from manifest
header_idxs = [i for i, c in enumerate(data.columns) if not c.endswith("file") and not c.startswith("histogram")]
header = list(data.columns[header_idxs])

# loop over rows
rows = []
for i in range(data.shape[0]):
    if i == 0:
        assert(data["compressor_name"][i] == opts.reference_row)
        continue

    # skip if not input files
    if not os.path.exists(data["spectra_file"][i]) or not os.path.exists(data["halo_finder_file"][i]):
        print("Input files do not exist! Skipping!")
        continue

    # get CSV data from manifest
    row = list(data.values[i][header_idxs])
    output_files = []

    # plot spectrum ratio
    output_files.append("spectrum_{}.png".format(i))
    if not os.path.exists(output_files[-1]):
        cmd = ["python", "plot_spectrum.py",
               "--input-file", data["spectra_file"][i],
               "--reference-file", data["spectra_file"][0],
               "--output-file", opts.output_file + "/" + output_files[-1],
               "--xlim", 0, 10,
               "--ylim", 0.99, 1.01,
               "--operation", "ratio"]
        _external_call(cmd, debug=opts.debug)

    # plot halo mass distribution
    output_files.append("halo_{}.png".format(i))
    if not os.path.exists(output_files[-1]):
        cmd = ["python", "plot_gio_distribution.py",
               "--input-file", data["halo_finder_file"][i], data["halo_finder_file"][0],
               "--output-file", opts.output_file + "/" + output_files[-1],
               "--xlim", 1e10, 1e15,
               "--ylim", 1e-1, 1e6,
               "--xlog",
               "--ylog",
               "--log-bins"]
        _external_call(cmd, debug=opts.debug)

    # plot halo mass distribution ratio
    output_files.append("halo_ratio_{}.png".format(i))
    if not os.path.exists(output_files[-1]):
        cmd = ["python", "plot_gio_distribution.py",
               "--input-file", data["halo_finder_file"][i],
               "--reference-file", data["halo_finder_file"][0],
               "--output-file", opts.output_file + "/" + output_files[-1],
               "--xlim", 1e10, 1e14,
               "--ylim", 0.5, 2.0,
               "--xlog",
               "--log-bins",
               "--operation", "ratio"]
        _external_call(cmd, debug=opts.debug)

    # plot halo mass distribution ratio
    output_files.append("halo_ratio_zoom_{}.png".format(i))
    if not os.path.exists(output_files[-1]):
        cmd = ["python", "plot_gio_distribution.py",
               "--input-file", data["halo_finder_file"][i],
               "--reference-file", data["halo_finder_file"][0],
               "--output-file", opts.output_file + "/" + output_files[-1],
               "--xlim", 1e10, 1e14,
               "--ylim", 0.9, 1.1,
               "--xlog",
               "--log-bins",
               "--operation", "ratio"]
        _external_call(cmd, debug=opts.debug)

    # plot halo mass distribution absolute value
    output_files.append("halo_abs_{}.png".format(i))
    if not os.path.exists(output_files[-1]):
        cmd = ["python", "plot_gio_distribution.py",
               "--input-file", data["halo_finder_file"][i],
               "--reference-file", data["halo_finder_file"][0],
               "--output-file", opts.output_file + "/" + output_files[-1],
               "--xlim", 1e10, 1e14,
               "--ylim", 1e-1, 1e3,
               "--xlog",
               "--ylog",
               "--log-bins",
               "--operation", "abs"]
        _external_call(cmd, debug=opts.debug)

    # plot halo mass distribution difference
    output_files.append("halo_diff_{}.png".format(i))
    if not os.path.exists(output_files[-1]):
        cmd = ["python", "plot_gio_distribution.py",
               "--input-file", data["halo_finder_file"][i],
               "--reference-file", data["halo_finder_file"][0],
               "--output-file", opts.output_file + "/" + output_files[-1],
               "--xlim", 1e10, 1e14,
               "--ylim", -1e3, 1e3,
               "--xlog",
               "--log-bins",
               "--operation", "diff"]
        _external_call(cmd, debug=opts.debug)

    # read metrics file
    metrics = pandas.read_csv(data["metric_file"][i])

    # add to header
    if i == 1:
        init_metrics_header = list(metrics.columns)
        header += init_metrics_header
        header += ["FILE_{}".format(i) for i, _ in enumerate(output_files)]
    else:
        assert(all(init_metrics_header == metrics.columns))

    # append rows
    for j in range(metrics.shape[0]):
        rows.append(row + list(metrics.values[j]) + output_files)
        assert(len(rows[-1]) == len(header))

# write CSV file
with open(opts.output_file + "/data.csv", "w") as fp:
    fp.write(",".join(header) + "\n")
    for row in rows:
        fp.write(",".join(map(str, row)) + "\n")
