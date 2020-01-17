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

parser.add_argument('inference')
parser.add_argument('originalFile')
parser.add_argument('comparison_files', type=parseList)
parser.add_argument('directories', type=parseList)

args = parser.parse_args()

#inference = sys.argv[1]
#originalFile = sys.argv[2]
#comparison_files = sys.argv[3]
#directories = sys.argv[4]




for i in range(comparison_files.len()):
    os.mkdir(args.directory[i])
    cmd = "python " + args.inference + " " + args.originalFile + " " + args.comparison_files[i] + " " + args.directory[i]
    print(cmd)
    os.system(cmd)