#! /usr/bin/env python
""" Writes the distribution of halos from an SQL query.
"""
import sys
sys.path.append(sys.path[0]+"/..")
import argparse
import numpy
from pat import gioSqlite as gio_sqlite

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
parser.add_argument("--input-file", required=True)
parser.add_argument("--output-file", required=True)
parser.add_argument("--sqlite-file", default="/projects/exasky/visio/genericio/frontend/GenericIOSQLite.so")
parser.add_argument("--query", default="select fof_halo_mass from __TABLE__ ORDER BY fof_halo_mass")
parser.add_argument("--xlabel", default="Halo Mass")
parser.add_argument("--ylabel", default="Counts")
parser.add_argument("--xlim", nargs="+", type=float, default=[])
parser.add_argument("--bins", type=float, default=20)
parser.add_argument("--log-bins", action="store_true")
opts = parser.parse_args()

# read data
data = numpy.array(load_sqlite_data(opts.input_file, opts.query, opts.sqlite_file))

# update ranges
x_min = min(x_min, data.min()) if not len(opts.xlim) > 0 else opts.xlim[0]
x_max = max(x_max, data.max()) if not len(opts.xlim) > 1 else opts.xlim[1]

# set binning and range of histograms
# can do uniform in linear space or logarithmic space
if opts.log_bins:
    bins = numpy.logspace(numpy.log10(x_min), numpy.log10(x_max), num=opts.bins)
    bins_range = None
else:
    bins = opts.bins
    bins_range = (x_min, x_max)

# create histogram
hist, bin_edges = numpy.histogram(data, bins=bins, range=bins_range)
hist = numpy.hstack([(0), numpy.repeat(hist, 2), (0)])
bin_edges = numpy.hstack([(bin_edges[0], bin_edges[0]),
                          numpy.repeat(bin_edges[1:-1], 2),
                          (bin_edges[-1], bin_edges[-1])])

# save results
delimiter = ","
results = numpy.column_stack([bin_edges, hist])
header = delimiter.join(map(str, [opts.xlabel, opts.ylabel]))
numpy.savetxt(opts.output_file, results, header=header, delimiter=delimiter)

print("Done!")
