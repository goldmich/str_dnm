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

def process_kong_markers():
	
	marker_set = None
	with open("./kong_et_al_2002_markers.txt", 'r') as open_kong_markers:
		header = open_kong_markers.readline()
		marker_set = set([x.split('\t')[1] for x in open_kong_markers.readlines()])
	
	marker_dict = {}
	with gzip.open("./stsInfo2_hg38.txt.gz", "rb") as open_sts_info:
		header = open_sts_info.readline()
		for l in open_sts_info.readlines():
			l_split = l.decode('ascii').split('\t')
			curr_intersect = set(l_split[7].split(',')) & marker_set
			if curr_intersect: # there should only ever be a single overlap
				if len(curr_intersect) > 1:
					pdb.set_trace()
				marker_dict[l_split[1]] = [curr_intersect.pop()]
	
	with gzip.open("./stsMap_hg38.txt.gz", "r") as open_sts_file:
		header = open_sts_file.readline()
		for l in open_sts_file.readlines():
			l_split = l.decode('ascii').split('\t')
			if l_split[3] in marker_dict:
				if not l_split[0] in all_chrs:
					continue
				marker_dict[l_split[3]] += l_split[:3]
	
	# i'm writing the markers in an unordered way ... i don't think it'll matter but i can change it later if i need to
	with open('./kong_et_al_2002_markers_hg38.txt', 'w') as open_outfile:
		open_outfile.write('#chrom\tstart\tend\tucsc_name\tkong_et_al_2002_name\n')
		for i in marker_dict:
			# a few markers do not map over and a few map over twice
			if len(marker_dict[i]) != 4:
				continue
			open_outfile.write('\t'.join(marker_dict[i][1:] + [i, marker_dict[i][0]]))
			open_outfile.write('\n')
	
	return None


if __name__ == '__main__':
	
	process_kong_markers()
