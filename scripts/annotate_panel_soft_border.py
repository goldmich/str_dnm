#!/usr/bin/env python

import pysam
import pdb
import re
import numpy as np
import os
import csv
import gzip
from scipy.spatial.distance import hamming

def load_hg38_ref(chr):
	outlines = ''
	ref_file = '/net/harris/vol1/data/hg38/%s.fa.gz' % chr
	with gzip.open(ref_file, 'r') as open_file:
		open_file.readline() # header
		outlines = ''.join([x.decode('ascii').rstrip('\n') for x in open_file.readlines()])
	return outlines
	

outlines = []
# infile = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/sites_annotated_reptiming_dnm_counts_ru.txt'
infile = snakemake.input[0]
with open(infile, 'r') as open_panel:
	curr_chr = ''
	curr_ref = ''
	curr_line = []
	last_end = -1
	for line in open_panel.readlines():
		line_split = line.rstrip('\n').split('\t')
		if line_split[0] != curr_chr:
			curr_chr = line_split[0]
			curr_ref = load_hg38_ref(curr_chr)
			last_end = -1
		curr_start = int(line_split[1]) - 1
		curr_end = curr_start + len(line_split[2]) * int(line_split[3])
		curr_line = line_split
		if curr_start - last_end <= 1:
			outlines[-1][-2] = "T" # this may change
			curr_line.append('T')
		else:
			curr_line.append('F')
		last_end = curr_end
		curr_template = (line_split[2] * 6)[:6]
		curr_line.append(
			min(
			hamming([x for x in curr_ref[curr_start - 6: curr_start].upper()],
					[x for x in curr_template]),
			hamming([x for x in curr_ref[curr_end: curr_end + 6].upper()],
					[x for x in curr_template])
			)
		)
		outlines.append(curr_line)

outfile = snakemake.output[0]
# # outfile = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/sites_annotated_reptiming_dnm_counts_ru_borders.txt'
with open(outfile, 'w') as open_outfile:
	open_outfile.write('\n'.join(['\t'.join([str(x) for x in y]) for y in outlines]) + '\n')


