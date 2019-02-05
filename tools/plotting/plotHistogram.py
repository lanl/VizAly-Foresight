#!/usr/bin/python

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

def splitIntoParts(number, parts):
    return np.linspace(0, number, parts+1)[1:]

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


def plotHistogram(y, plotName, maxVal, liveDisplay):
	# Create the x values
	numVals = len(y)
	x = splitIntoParts(maxVal, numVals)

	# What to draw
	plt.plot(x,y)

	plt.title(plotName)
	plt.ylabel("Frequency")

	# Display
	plt.savefig(plotName)

	if (liveDisplay):
		plt.show()


def main(argv):
	# check if we have arguments
	if (len(sys.argv) < 3):
		print ("Arguments needed")
		sys.exit(0)

	# read csv file
	rows = readCSV(sys.argv[1])

	# draw
	plotHistogram(rows, sys.argv[2], float(sys.argv[3]), True)


if __name__ == "__main__":
    main(sys.argv)

# Usage
# python plotHistogram.py ../../build/BigCrunch_y_absolute_error_hist.csv BigCrunch_abs_err_x 0.124985