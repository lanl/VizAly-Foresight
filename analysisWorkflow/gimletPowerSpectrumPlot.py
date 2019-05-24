#!/usr/bin/python

import sys, json, os, csv
import matplotlib.pyplot as plt


def extractValue(filename, colpos):
	col = []
	with open(filename) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=' ')
		for row in readCSV:
			col.append( float(row[colpos]) )

	return col



def plotGraph(x, x_label, y_label, title, list_of_tuples):
	fig = plt.figure()
	ax = plt.subplot(1,1,1)

	plt.title(title)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.grid(True)

	for item in list_of_tuples:
		ax.plot(x, item[0], label=item[1], marker=item[2])
		ax.legend()

	fig = plt.gcf()
	plt.show()
	fig.savefig(title + '.png', dpi=100)



if __name__ == "__main__":
	if len(sys.argv) < 2:
		print ("Json file and title needed; e.g. python gimletPowerSpectrumPlot.py gimletInput.json")
		exit()

	print sys.argv[1]

	# Read input json file
	with open(sys.argv[1], "r") as read_file:
		json_data = json.load(read_file)


	for ana in json_data["analysis"]:
		plot_title = ana["title"]
		print plot_title

		to_plot = []  # all the items to plot

		k_list = []
		orig_pk = []
		for file in ana["files"]:
			if (file["name"]=="orig"):
				k_list  = extractValue(file["path"], 2)
				orig_pk = extractValue(file["path"], 3)
			else:
				temp_pk = extractValue(file["path"], 3)
				pk_ratio = [i / j for i, j in zip(temp_pk, orig_pk)]
				this_tuple = (pk_ratio, file["name"], ".") #array, name, marker
				to_plot.append(this_tuple)

		plotGraph(k_list, 'k', 'pk', plot_title, to_plot)

"""
python gimletPowerSpectrumPlot.py gimletInput.json
"""