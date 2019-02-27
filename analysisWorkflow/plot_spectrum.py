#! /usr/bin/env python

import argparse
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
    else:
        raise NotImplementedError("Do not understand operation {}!".format(operation))

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument("--input-files", nargs="+", required=True)
parser.add_argument("--reference-file")
parser.add_argument("--output-file")
parser.add_argument("--operation", choices=["none", "ratio"], default="none")
parser.add_argument("--ylim", nargs="+", type=float, default=None)
parser.add_argument("--xlim", nargs="+", type=float, default=None)
parser.add_argument("--ylabel", default="P(k)")
parser.add_argument("--xlabel", default="k")
parser.add_argument("--ylog", action="store_true")
parser.add_argument("--xlog", action="store_true")
opts = parser.parse_args()

# load reference data
if opts.reference_file:
    ref = numpy.loadtxt(opts.reference_file, usecols=[0, 1])
    k_ref = ref[:, 0]
    pk_ref = ref[:, 1]

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
for input_file in opts.input_files:

    # read data
    data = numpy.loadtxt(input_file, usecols=[0, 1])
    k = data[:, 0]
    pk = data[:, 1]
    assert(numpy.all(k == k_ref))

    # only keep values in plot canvas
    idxs = []
    if opts.xlim:
        if len(opts.xlim) == 1:
            idxs += numpy.where(k > opts.xlim[0])
        else:
            idxs += numpy.where((k > opts.xlim[0]) & (k < opts.xlim[1]))
    if opts.ylim:
        if len(opts.ylim) == 1:
            idxs += numpy.where(pk > opts.ylim[0])
        else:
            idxs += numpy.where((pk > opts.ylim[0]) & (pk < opts.ylim[1]))
    idxs = sum([list(arr) for arr in idxs], [])
    idxs = numpy.unique(idxs)
    k = k[idxs]
    pk = pk[idxs]

    if opts.reference_file:
        pk_tmp = pk_ref[idxs]
        pk = operate(pk, pk_tmp, opts.operation)

    # plot
    plt_func(k, pk, label=os.path.basename(input_file))

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

if not opts.output_file:
    plt.show()
else:
    plt.savefig(opts.output_file)
