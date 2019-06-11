#!/usr/bin/python

import matplotlib
import itertools
matplotlib.use('agg')

import matplotlib.pyplot as plt


# From: https://stackoverflow.com/questions/29461608/matplotlib-fixing-x-axis-scale-and-autoscale-y-axis
def autoscale_y(ax,margin=0.1):
    import numpy as np

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo,hi = ax.get_xlim()
        y_displayed = yd[((xd>lo) & (xd<hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed)-margin*h
        top = np.max(y_displayed)+margin*h
        return bot,top

    lines = ax.get_lines()
    bot,top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot: bot = new_bot
        if new_top > top: top = new_top

    ax.set_ylim(bot,top)



def plotGraph(x, x_label, y_label, title, list_of_tuples):
	fig = plt.figure()
	ax = plt.subplot(1,1,1)

	plt.title(title)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.grid(True)

	marker = itertools.cycle(('o', 's', '*', 'P', 'D', 'X', 'v', 'H'))

	for item in list_of_tuples:
		ax.semilogx(x, item[0], label=item[1], marker=marker.next())
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.00), shadow=True, ncol=2)

	plt.xlim(0, 12)

	

	fig = plt.gcf()
	autoscale_y(ax)
	plt.show()
	fig.savefig(title + '.png', dpi=100)
