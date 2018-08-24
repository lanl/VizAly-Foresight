#!/usr/bin/python
import sys
import csv
import json
import random
import itertools
import matplotlib.pyplot as plt



def drawScatter(oufputFileName, compressorsToPlot, metricsToPlot, rows, fields):

	#Find rows to plot
	dataRowLocation = []
	for c in compressorsToPlot:
		index = 0
		for r in rows:
			if (r[0] == c):
				dataRowLocation.append(index)
			index = index + 1


	# find value location
	metricIndexToPlot = []
	for m in metricsToPlot:
		
		index = 0
		for f in fields:
			if (f.lstrip() == m):
				metricIndexToPlot.append(index)
			index = index + 1


	# Draw plot
	color  = itertools.cycle( ('b', 'g', 'c', 'y', 'm', 'r', 'k') )
	marker = itertools.cycle( (',', 'x', '+', 'v', 'o', '*', '+') ) 

	lines = []
	for i in range( len(dataRowLocation)):
		x = plt.scatter( float(rows[ dataRowLocation[i] ][ metricIndexToPlot[0] ]), float(rows[ dataRowLocation[i] ][ metricIndexToPlot[1] ]), marker=marker.next(), color=color.next())
		lines.append(x)
		index = index + 1

	plt.xlabel(metricsToPlot[0])
	plt.ylabel(metricsToPlot[1])

	plt.legend(lines,
           compressorsToPlot,
           scatterpoints=1,
           loc='best',
           ncol=2,
           fontsize=10)

	plt.savefig(oufputFileName)
	plt.show()



def drawBar(oufputFileName, compressorsToPlot, metricToPlot, rows, fields):	

	# Find rows to plot
	dataRowLocation = []
	for c in compressorsToPlot:
		
		index = 0
		for r in rows:
			if (r[0] == c):
				dataRowLocation.append(index)
			index = index + 1

	# find value location
	index = 0
	metricIndexToPlot = -1
	for f in fields:
		if (f.lstrip() == metricToPlot):
			metricIndexToPlot = index

		index = index + 1


	# Fill in data
	y = []
	for i in dataRowLocation:
		y.append( float(rows[i][metricIndexToPlot]) )

	x = range( len(y) )


	# naming the axes
	plt.xlabel('Compressor_field')
	plt.ylabel(metricToPlot)

	# plotting a bar chart
	plt.bar(x, y, tick_label = compressorsToPlot, width = 0.8, color = ['red', 'green'])


	plt.savefig(oufputFileName)
	plt.show()


# Read a csv file
def readCSV(inputFileName, fields, rows):
	with open(inputFileName, 'r') as csvfile:
		csvreader = csv.reader(csvfile)

		# read header
		fields = csvreader.next()

		# read data
		for row in csvreader:
			rows.append(row)

	return fields, row


def main(argv):

	# check if we have arguments
	if (len(sys.argv) < 1):
		print ("Arguments needed")
		sys.exit(0)

	# Open json file
	with open(sys.argv[1]) as jsonFile:
		jsonData = json.load(jsonFile)

	# parse command line arguments
	inputFileName = jsonData["inputFileName"]


	#initialize data storage
	fields = []	# field names
	rows = []	# data

	# read file
	fields, row = readCSV(inputFileName, fields, rows)


	compressorsToPlot = []
	for c in jsonData["x-axis"]:
		compressorsToPlot.append(c)

	#metricToPlot = "Compression Ratio"
	#metricsToPlot = ["Compression Throughput", "Compression Ratio"]


	if (jsonData["plot-type"] == "scatter"):
		metricsToPlot = jsonData["metrics"]
		drawScatter(jsonData["outputFileName"], compressorsToPlot, metricsToPlot, rows, fields)

	if (jsonData["plot-type"] == "barchart"):
		metricToPlot = jsonData["metrics"]
		drawBar(jsonData["outputFileName"], compressorsToPlot, metricToPlot, rows, fields)
	
	
	#
	
	

if __name__ == "__main__":
    main(sys.argv)

# Requirements:
#  - python 2

# Setup:
#  - sudo apt install python-pip
#  - sudo apt-get install python-setuptools
#  - pip install matplotlib

# sudo apt-get install python-tk

# Run:
# python plotting.py metrics
