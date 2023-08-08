#!/usr/bin/env python

import csv
import glob
import re
import gzip
import pdb
import argparse
import os
import math
import itertools
import random
import numpy as np

autosomes = ["chr" + str(x) for x in range(1, 23)]
all_chrs = autosomes + ['chrX']

def annotate_file_kong_markers(str_file, kong_markers, outfile_name):
	
	kong_dict = {}
	with open(kong_markers, 'r') as open_kong_file:
		header = open_kong_file.readline()
		for l in open_kong_file:
			l_split = l.rstrip('\n').split('\t')
			if l_split[0] not in kong_dict:
				kong_dict[l_split[0]] = []
			kong_dict[l_split[0]].append((int(l_split[1]), int(l_split[2])))
	
	with open(outfile_name, 'w') as open_outfile:
		with open(str_file, 'r') as open_panel_file:
			counter = 0
			outlines = ''
			for l in open_panel_file:
				l_split = l.rstrip('\n').split('\t')
				l_split[1] = int(l_split[1])
				if [x for x in kong_dict[l_split[0]] if x[0] <= l_split[1] and x[1] >= l_split[1]]:
					outlines += '\n' + l.rstrip('\n') + '\tT'
				else:
					outlines += '\n' + l.rstrip('\n') + '\tF'
				counter += 1
				if counter > 1000:
					counter = 0
					open_outfile.write(outlines)
					outlines = ''
			# remaining outlines
			open_outfile.write(outlines)

	return None

if __name__ == '__main__':
	
	annotate_file_kong_markers(snakemake.input[0], snakemake.input[1], snakemake.output[0])
	
