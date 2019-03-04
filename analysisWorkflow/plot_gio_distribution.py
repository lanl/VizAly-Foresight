#! /usr/bin/env python
""" Plots the distribution of halos from an SQL query.
"""

import argparse
import gioSqlite as gio_sqlite
import matplotlib.pyplot as plt
import numpy
import os

def operate(y, ref, operation):
    """ Suite of math operations for plot.
    """
    if operation == "none":
        return y
    elif operation == "ratio":
        return y / ref
    elif operation == "abs":
        return numpy.abs(y - ref)
    elif operaiton == "diff":
        return y - ref
    else:
        raise NotImplementedError("Do not understand operation {}!".format(operation))

def load_sqlite_data(path, query, sqlite_file):
    """ Loads data using SQLite query.
    """

    # load file
    print("Reading {}...".format(path))
    query_mgr = gio_sqlite.GioSqlite3()
    query_mgr.loadGIOSqlite(sqlite_file)

    # load data
    i = 0
    table_name = "foo_{}".format(i)
    query_mgr.createTable(table_name, (path))

    # execute query
    query = query.replace("__TABLE__", table_name)
    result = query_mgr.runQueryOutputList(query)

    # typecast
    result = numpy.array(result).flatten()
    assert(len(result.shape) == 1)

    return result

# parse command line
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input-files", nargs="+", required=True)
parser.add_argument("--reference-file")
parser.add_argument("--output-file")
parser.add_argument("--sqlite-file", default="/home/cmbiwer/src/VizAly-Vis_IO/genericio/frontend/GenericIOSQLite.so")
parser.add_argument("--query", default="select fof_halo_mass from __TABLE__ ORDER BY fof_halo_mass DESC LIMIT 50")
parser.add_argument("--operation", choices=["none", "ratio", "abs", "diff"], default="none")
parser.add_argument("--ylim", nargs="+", type=float, default=[])
parser.add_argument("--xlim", nargs="+", type=float, default=[])
parser.add_argument("--ylabel", default="Counts")
parser.add_argument("--xlabel", default="Halo Mass")
parser.add_argument("--labels", nargs="+", default=None)
parser.add_argument("--ylog", action="store_true")
parser.add_argument("--xlog", action="store_true")
parser.add_argument("--log-bins", action="store_true")
parser.add_argument("--bins", type=int, default=25)
opts = parser.parse_args()

# sanity check labels
if opts.labels:
    assert(len(opts.input_files) == len(opts.labels))

# load reference data
if opts.reference_file:
    ref = load_sqlite_data(opts.reference_file, opts.query, opts.sqlite_file)

# determine plot function
if opts.ylog and opts.xlog:
    plt_func = plt.loglog
elif opts.ylog:
    plt_func = plt.semilogy
elif opts.xlog:
    plt_func = plt.semilogx
else:
    plt_func = plt.plot

# load input files data
x_min = numpy.inf
x_max = 0.0
data = {}
for i, input_file in enumerate(opts.input_files):

    # read data
    tmp = load_sqlite_data(input_file, opts.query, opts.sqlite_file)

    # store results
    label = opts.labels[i] if opts.labels else os.path.basename(input_file).strip("__")
    data[label] = numpy.array(tmp)

    # update range
    x_min = min(x_min, data[label].min()) if not len(opts.xlim) > 0 else opts.xlim[0]
    x_max = max(x_max, data[label].max()) if not len(opts.xlim) > 1 else opts.xlim[1]

# create histograms on equal x-axis range and plot
keys = list(data.keys())
keys.sort()
for i, key in enumerate(keys):

    # set binning and range of histograms
    if opts.log_bins and i == 0:
        bins = numpy.logspace(numpy.log10(x_min), numpy.log10(x_max), num=opts.bins)
        bins_range = None
    elif i == 0:
        bins = opts.bins
        bins_range = (x_min, x_max)

    # create histogram
    hist, bin_edges = numpy.histogram(data[key], bins=bins, range=bins_range)
    hist = numpy.hstack([(0), numpy.repeat(hist, 2), (0)])
    bin_edges = numpy.hstack([(bin_edges[0], bin_edges[0]),
                              numpy.repeat(bin_edges[1:-1], 2),
                              (bin_edges[-1], bin_edges[-1])])

    # operate on histogram with reference file
    if opts.reference_file:
        if i == 0:
            ref, _ = numpy.histogram(ref, bins=bins, range=bins_range)
            ref = numpy.hstack([(0), numpy.repeat(ref, 2), (0)])
        hist = operate(hist, ref, opts.operation)

    # plot
    plt_func(bin_edges, hist, label=key)

# format plot
if opts.xlim:
    xlim = opts.xlim[0] if len(opts.xlim) == 1 else opts.xlim
    plt.xlim(xlim)
if opts.ylim:
    ylim = opts.ylim[0] if len(opts.ylim) == 1 else opts.ylim
    plt.ylim(ylim)
plt.xlabel(opts.xlabel)
plt.ylabel(opts.ylabel)
plt.grid()
plt.legend()

# display
if not opts.output_file:
    plt.show()
else:
    plt.savefig(opts.output_file)

