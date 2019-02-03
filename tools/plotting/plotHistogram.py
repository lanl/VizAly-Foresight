#!/usr/bin/python

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

# Read a csv file
def readCSV(inputFileName):
	rows = []
	with open(inputFileName, 'r') as csvfile:
		csvreader = csv.reader(csvfile)

		# read data
		for row in csvreader:
			rows.append(row)

	# Remove last item which is empty
	row = row[:-1]

	# Convert to float
	myList = [ float(x) for x in row ]
	return myList


def plotHistogram(val, plotName, oufputFileName):
	x = range( len(val) )

	plt.bar(x, val, width = 1.0, align='edge')

	plt.title(plotName)
	plt.xlabel("Absolute Error")
	plt.ylabel("Frequency")
	plt.xticks([])

	liveDisplay = false
	if (liveDisplay):
		plt.show()

	plt.savefig(oufputFileName)


def main(argv):
	# check if we have arguments
	if (len(sys.argv) < 3):
		print ("Arguments needed")
		sys.exit(0)

	# read csv file
	rows = readCSV(sys.argv[1])

	# draw
	plotHistogram(rows, sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main(sys.argv)