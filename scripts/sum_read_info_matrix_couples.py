#!/usr/bin/env python

import random
import pysam
import pdb
import re
import numpy as np
import os

# matrix_filename_list = snakemake.input[0:]
matrix_filename_list = ['/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/12220.22.parents_callable_dnms_matrix.txt',
						'/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/12220.22.parents_callable_dnms_matrix.txt']
family_matrix = None
curr_couple = None
header = None

for curr_file in matrix_filename_list:
	if not curr_couple:
		curr_matrix = np.loadtxt(curr_file, dtype = str)
		header = np.reshape(curr_matrix[0,], (1, 7))
		curr_couple = (curr_matrix[1,0], curr_matrix[1,1])
		curr_matrix = curr_matrix[1:,2:].astype(int)
		family_matrix = curr_matrix
	else:
		curr_matrix = np.loadtxt(curr_file, skiprows = 1, usecols = (2, 3, 4, 5, 6), dtype = int)
		family_matrix[0:,4] += curr_matrix[0:, 4]

couple_name_matrix = np.empty((family_matrix.shape[0], 2), dtype = '<U11')
couple_name_matrix[:, 0] = [curr_couple[0] for _ in range(couple_name_matrix.shape[0])]
couple_name_matrix[:, 1] = [curr_couple[1] for _ in range(couple_name_matrix.shape[0])]
family_matrix_cat = np.concatenate((couple_name_matrix, family_matrix), axis = 1)
matrix_to_output = np.concatenate((header, family_matrix_cat), axis = 0)

output_filename = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/12220.test.parents_callable_dnms_matrix.txt'
# output_filename = snakemake.output[0]
with open(output_filename, 'w') as open_outfile:
	for l in matrix_to_output:
		open_outfile.write('\t'.join(l) + '\n')

