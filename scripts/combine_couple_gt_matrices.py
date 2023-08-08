#!/usr/bin/env python

import pdb
import re
import numpy as np
import os

# curr_phase = 'phase3_1'

couples_dict = {}
header = ''
for infile in snakemake.input[:23]:
# for infile in ["../tmp/phase3_1_1.indiv_repeat_content.txt", "../tmp/phase3_1_2.indiv_repeat_content.txt"]:
# 	print(infile)
	with open(infile, 'r') as open_matrix_file:
		header = open_matrix_file.readline()
		for l in open_matrix_file.readlines():
			l_split = l.rstrip('\n').split('\t')
			couple = tuple(l_split[:2])
			if couple not in couples_dict:
				couples_dict[couple] = np.zeros((4, 3), dtype = int)
			couples_dict[couple] = couples_dict[couple] + np.reshape([int(x) for x in l_split[2:]], (4, 3))

with open(snakemake.output[0], 'w') as open_outfile:
	open_outfile.write(header)
	for s1, s2 in couples_dict:
		open_outfile.write(s1 + '\t' + s2 + '\t')
		open_outfile.write('\t'.join([str(x) for x in np.reshape(couples_dict[(s1, s2)][0:,], -1)]) + '\n')
