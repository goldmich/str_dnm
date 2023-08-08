#!/usr/bin/env python

import pdb
import re
import numpy as np
import os

# curr_phase = 'phase3_1'

sample_dict = {}
header = ''
for infile in snakemake.input[:23]:
# for infile in ["../tmp/phase3_1_1.indiv_repeat_content.txt", "../tmp/phase3_1_2.indiv_repeat_content.txt"]:
# 	print(infile)
	with open(infile, 'r') as open_matrix_file:
		header = open_matrix_file.readline()
		for l in open_matrix_file.readlines():
			l_split = l.rstrip('\n').split('\t')
			if l_split[0] not in sample_dict:
				sample_dict[l_split[0]] = np.zeros((4, 2), dtype = int)
			sample_dict[l_split[0]] = sample_dict[l_split[0]] + np.reshape([int(x) for x in l_split[1:]], (4, 2))

with open(snakemake.output[0], 'w') as open_outfile:
	open_outfile.write(header)
	for s in sample_dict:
		open_outfile.write(s + '\t')
		open_outfile.write('\t'.join([str(x) for x in np.reshape(sample_dict[s][0:,], -1)]) + '\n')
