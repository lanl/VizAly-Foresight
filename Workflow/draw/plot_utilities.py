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



def plotScatterGraph(opts, x_data, y_data, label):
    # Type of plot
    if opts.ylog and opts.xlog:
        plt_func = plt.loglog
    elif opts.ylog:
        plt_func = plt.semilogy
    elif opts.xlog:
        plt_func = plt.semilogx
    else:
        plt_func = plt.plot


    plt_func(x_data, y_data, label)


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

    # display + print
    plt.show()
    plt.savefig(opts.output_file)



def plotScatterGraph(x, x_label, y_label, title, path, x_range, list_of_tuples):
    fig = plt.figure()
    ax = plt.subplot(1,1,1)

    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(True)

    marker = itertools.cycle(('o', 's', '*', 'P', 'D', 'X', 'v', 'H', '2', '+'))

    for item in list_of_tuples:
        ax.semilogx(x, item[0], label=item[1], marker=next(marker))
        ax.legend(bbox_to_anchor=[1.05,1],loc='upper left',borderaxespad=0)

    plt.xlim(x_range[0], x_range[1])

    fig = plt.gcf()
    autoscale_y(ax)
    plt.show()
    fig.savefig(path +'/' + title + '.png', dpi=100, bbox_inches="tight")


def create_img(text, filename):
	fig, ax = plt.subplots() 
   
	fig.text(0.2, 0.4, text, fontsize=30, color="black") 
	ax.set(xlim = (0, 8), ylim = (0, 8)) 

	_filename = filename.replace(' ','_')
	plt.savefig(_filename + ".png")


def generic_plot(x, y_orig, x_label, y_label, name, y_tuple_labels=[], plot_type="linear", axis_type="simple"):
    """
        axis_type options: linear, x_log, y_log, log
        plot_type:  simple (one y value)
                    two (plot both y values)
                    diff (multiple y values: y_0 - y_1)
                    ratio (multiple y values: y_0/y_1)
    """

    # Extract Data
    num_files = len(y_orig)


    # Process data if needed
    y = []
    if num_files > 1:
        if plot_type == "linear":
            y = y_orig

        elif plot_type == "two":
            for i in range(num_files):
                y.append( y_orig[i] )

        elif plot_type == "diff":
            
            for i in range(num_files-1):
                temp_y = []
                for j in range( len(y_orig[0]) ):
                    temp_y.append(  y_orig[i+1][j] - y_orig[0][j] )
                y.append( temp_y )

                
        elif plot_type == "ratio":
            
            for i in range(num_files-1):
                temp_y = []
                for j in range( len(y_orig[0]) ):
                    if ( y_orig[i+1][j] == 0):
                        temp_y.append( (y_orig[i+1][j]+1) / (y_orig[0][j]+1) )
                    else:
                        temp_y.append( y_orig[i+1][j] / y_orig[0][j] )
                y.append( temp_y )
    else:
        y = y_orig


    # Create layout
    color_list = ["blue", "red", "orange", "green",  "black"]

    fig = plt.figure()
    ax = plt.subplot(1,1,1)

    plt.title(name)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(True)

    
    for i in range(len(y)):
        line_label=""
        try:
            line_label = y_tuple_labels[i]
        except:
            line_label=""
        plt.plot(x, y[i], '-', color=color_list[i], label=line_label, marker='')

    if axis_type == "x_log":
        plt.xscale('log')
    elif axis_type == "y_log":
        plt.yscale('log')
    elif axis_type == "log":
        plt.xscale('symlog')
        plt.yscale('symlog')
    
    plt.legend()
    filename = name.replace(' ','_')
    plt.savefig(filename + ".png")



def read_files(filenames, seperator=',', x_col=0, y_col=1, limit_type="None", limit_val=1):
    """
        limit_type: 
            None: read all
            x_less: x less than limit_value
    """

    num_files = len(filenames)

    x = extract_csv_col(filenames[0], seperator, x_col)
    
    if x == None:
        x = []
        y = []
        return x, y
    

    y_orig = []
    for i in range(num_files):
        y_orig.append( extract_csv_col(filenames[i], seperator, y_col) )


    if limit_type == "x_less":
        count = 0
        for i in x:
            if i < limit_val:
                count = count + 1

        new_x = x[:count]
        new_y = []
        for i in range(num_files):
            new_y.append(y_orig[i][:count])

        return new_x, new_y
    else:
        return x, y_orig


def create_plot(filenames, x_label, y_label, name, seperator=',', x_col=0, y_col=1, limit_type="None", limit_val=1, y_tuple_labels=[], plot_type="linear", axis_type="simple"):
	x, y = read_files(filenames, seperator, x_col, y_col, limit_type, limit_val)
	if x == [] or y == []:
		print("Could not create", name)      
		create_img("Can't create\n"+name, name)
	else:
		generic_plot(x, y, x_label, y_label, name, y_tuple_labels, plot_type, axis_type)

