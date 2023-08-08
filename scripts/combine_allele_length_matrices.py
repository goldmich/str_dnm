#!/usr/bin/env python

import pdb
import re
import numpy as np
import os

# curr_chr = 'chr1'

allele_lengths_dict = {}
header = ''
for infile in snakemake.input[:6]:
# for infile in ["../tmp/phase3_1_10.allele_lengths.txt", "../tmp/phase3_2_10.allele_lengths.txt"]:
# 	print(infile)
	with open(infile, 'r') as open_matrix_file:
		header = open_matrix_file.readline()
		for l in open_matrix_file.readlines():
			l_split = l.rstrip('\n').split('\t')
			curr_locus = (l_split[0], int(l_split[1]))
			if curr_locus not in allele_lengths_dict:
				allele_lengths_dict[curr_locus] = np.zeros((100), dtype = int)
			allele_lengths_dict[curr_locus] = np.logical_or(allele_lengths_dict[curr_locus], [int(x) for x in l_split[2:]])

# with open('test_combine_allele_lengths.txt', 'w') as open_outfile:
with open(snakemake.output[0], 'w') as open_outfile:
	open_outfile.write(header)
	for chr, pos in sorted(allele_lengths_dict.keys(), key=lambda x:x[1]):
		open_outfile.write(chr + '\t' + str(pos) + '\t')
		open_outfile.write('\t'.join([str(int(x)) for x in allele_lengths_dict[(chr, pos)]]) + '\n')
