#! /usr/bin/env python
""" Plots distributions of hashed GenericIO files.
"""

import argparse
import gioSqlite as gio_sqlite
import matplotlib.pyplot as plt
import numpy
import os

def load_sqlite_data(path, query, sqlite_file):

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
    result = numpy.array(result)

    return result

def snap(x, d=0.25):
#    d = 1.0 / d
#    return numpy.ceil(d * x) / d
    return numpy.ceil(x / d) * d

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=int, default=0)
parser.add_argument("--parameters", nargs="+", default=["fof_halo_center_x", "fof_halo_center_y", "fof_halo_center_z", "fof_halo_mass"])
parser.add_argument("--match-parameters", nargs="+", default=["fof_halo_center_x"])
parser.add_argument("--match-parameters-snaps", nargs="+", type=float, default=[-1])
parser.add_argument("--order-by-parameter", default="fof_halo_mass")
parser.add_argument("--input-file")
parser.add_argument("--reference-file")
parser.add_argument("--output-file-format")
parser.add_argument("--sqlite-file", default="/home/cmbiwer/src/VizAly-Vis_IO/genericio/frontend/GenericIOSQLite.so")
parser.add_argument("--limit", type=float, default=None)
parser.add_argument("--bins", type=int, default=50)
opts = parser.parse_args()

# add parameters
for p in opts.match_parameters:
    if p not in opts.parameters:
        opts.parameters.append(p)

# create query
#query = "select {} from __TABLE__ ORDER BY {} DESC".format(", ".join(opts.parameters), opts.order_by_parameter)
query = "select {} from __TABLE__".format(", ".join(opts.parameters))
if opts.limit:
    query += " LIMIT {}".format(opts.limit)

# set seed
numpy.random.seed(opts.seed)

# count scalar fields to query and match
n_parameters = len(opts.parameters)
assert(len(opts.match_parameters) == len(opts.match_parameters_snaps))

# determine index of match parameters
match_parameters_idxs = []
for p in opts.match_parameters:
    match_parameters_idxs.append(opts.parameters.index(p))

# set data to compare
if opts.input_file and opts.reference_file:
    data_1 = load_sqlite_data(opts.reference_file, query, opts.sqlite_file)
    data_2 = load_sqlite_data(opts.input_file, query, opts.sqlite_file)
    data_1 = numpy.transpose(data_1)
    data_2 = numpy.transpose(data_2)
else:
    raise KeyError("Must use both --input-file and --reference-file!")

# print statment
print("Read {} particles from {}...".format(data_1.shape[1], opts.input_file))
print("Read {} particles from {}...".format(data_2.shape[1], opts.reference_file))

# snap
data_1_snap = []
data_2_snap = []
for i, d in zip(match_parameters_idxs, opts.match_parameters_snaps):
    if d == -1:
        data_1_snap.append(data_1[i, :])
        data_2_snap.append(data_2[i, :])
    else:
        data_1_snap.append(snap(data_1[i, :], d))
        data_2_snap.append(snap(data_2[i, :], d))
    print(data_1[i, :])
    print(data_1_snap[-1])
data_1_snap = numpy.vstack(data_1_snap)
data_2_snap = numpy.vstack(data_2_snap)

# transpose for loop
data_1_snap = numpy.transpose(data_1_snap)
data_2_snap = numpy.transpose(data_2_snap)
data_1_snap[data_1_snap == -0.0] = 0.0
data_2_snap[data_2_snap == -0.0] = 0.0

print
print(data_1)
print(data_2)
print(data_1_snap)
print(data_2_snap)

# hash
n_lost = 0
diffs = []
lost_idxs = []
for i, d_1 in enumerate(data_1_snap):

    # find matches
    idx = numpy.flatnonzero((d_1 == data_2_snap).all(1))

    # ensure unqiue matches
    # record number of particles without match
    if idx.size > 1:
        #print("Too many FOUND!", i, idx.size)
        n_lost += 1
        continue
#        raise KeyError("Mapping to more than one particle!")
    elif idx.size == 0:
        n_lost += 1
        lost_idxs.append(i)
        #print("Not FOUND!", i, idx)
        continue
    else:
        idx = idx[0]

    # get hashed values
    d_2 = data_2_snap[idx]

    # get data values
    r_1 = data_1[:, i]
    r_2 = data_2[:, idx]

    # compute differences
    diff = r_2 - r_1
    diff_hash = d_2 - d_1

    # append result
    diffs.append(diff)

# typecast and count
diffs = numpy.vstack(diffs)
n_found = diffs.shape[0]
n_particles = data_1.shape[1]

# print statement
print("Number of matched particles: {}".format(n_found))
print("Number of missed particles: {}".format(n_lost))
print("Indices missed were: {}".format(lost_idxs))
#for lost_idx in lost_idxs:
#    print()
#    print(data_1[:, lost_idx])
#    print(data_1_snap[lost_idx, :])
assert(n_found + n_lost == n_particles)

# plot
for i in range(n_parameters):
    plt.hist(diffs[:, i], bins=opts.bins)
    plt.xlabel("Difference {}".format(opts.parameters[i]))
    if opts.output_file_format:
        plt.savefig(opts.output_file_format.format(i))
    else:
        plt.show()
    plt.close()
