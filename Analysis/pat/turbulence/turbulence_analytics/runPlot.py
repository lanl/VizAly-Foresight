import os
import sys
import re
import argparse


def parseList(s):
	# s should be in the format ['aa','bb','cc']
	try:
		s1 = s.replace('[','')
		s2 = s1.replace(']','')
		s3 = s2.replace('\'','')

		l = list(map(str, s3.split(',')))
		return l

	except:
		raise argparse.ArgumentTypeError("Something is wrong with the list")



# Gather inputs
parser = argparse.ArgumentParser()

parser.add_argument('plot')
parser.add_argument('data_directories', type=parseList)
parser.add_argument('ourput_prefixes', type=parseList)
parser.add_argument('plot_dir')

args = parser.parse_args()

#plot = sys.argv[1]
#directories = sys.argv[2]
#ourput_prefixes = sys.argv[3]
#plot_dir = sys.argv[4]


os.mkdir(plot_dir)

for i in range(directories.len()):
    cmd = "python " + args.plot + " " + args.data_directories[i] + " " + args.ourput_prefixes[i] + " " + args.plot_dir
    print(cmd)
    os.system(cmd)