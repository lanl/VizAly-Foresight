#!/usr/bin/python

import sys, json, os, csv
import matplotlib.pyplot as plt

filename_orig = "/home/pascal/Desktop/gimletPlotting/orig_rhob_ps3d.txt"
filename_sz = "/home/pascal/Desktop/gimletPlotting/sz_rhob_ps3d.txt"
title_name = "ratio_rhob_ps3d"


def extractValue(filename, colpos):
	col = []
	with open(filename) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=' ')
		for row in readCSV:
			col.append( float(row[colpos]) )

	return col


def plotGraph(x, x_label, y_label, title, list_of_tuples):
	fig = plt.figure()
	ax = plt.subplot(111)

	plt.title(title)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.grid(True)

	for item in list_of_tuples:
		ax.plot(x, item[0], label=item[1], marker=item[2])
		ax.legend()
	
	fig = plt.plot(figsize=(50,50))
	plt.show()



if __name__ == "__main__":
	k_list  = extractValue(filename_orig, 2)
	orig_pk = extractValue(filename_orig, 3)
	sz_pk   = extractValue(filename_sz, 3)

	sz_pk_ratio = [i / j for i, j in zip(sz_pk, orig_pk)] 

	toplot = []

	sz = (sz_pk_ratio, "SZ", ".") #array, name, marker
	toplot.append(sz)

	plotGraph(k_list, 'k', 'pk', 'rhob', toplot)
