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

def annotate_file_with_rep_timing_and_cpg(str_file, allele_lengths_file, outfile_name):
	
	
	allele_lengths = {}
	with open(allele_lengths_file, 'r') as open_allele_lengths:
		header = open_allele_lengths.readline()
		for l in open_allele_lengths:
			l_split = l.rstrip('\n').split('\t')
			if l_split[0] not in allele_lengths:
				allele_lengths[l_split[0]] = {}
			allele_lengths[l_split[0]][int(l_split[1])] = [int(x) for x in l_split[2:]]
	
	with open(outfile_name, 'w') as open_outfile:
		with open(str_file, 'r') as open_dnm_file:
			# header
			header = open_dnm_file.readline().rstrip('\n')
			open_outfile.write(header + ',rsdnm')
			counter = 0
			outlines = ''
			for l in open_dnm_file:
				l_split = l.rstrip('\n').split(',')
				# only autosomes
				if l_split[0] not in allele_lengths:
					outlines += '\n' + l.rstrip('\n') + ',-1'
					continue
				# python is 0-indexed, of course, but i have kinda 1-indexed this allele matrix (to allow for hypothetically alleles of length 0)
				outlines += '\n' + l.rstrip('\n') + ',' + str(allele_lengths[l_split[0]][int(l_split[1])][int(l_split[8])])
				counter += 1
				if counter > 1000:
					counter = 0
					open_outfile.write(outlines)
					outlines = ''
			# remaining outlines
			open_outfile.write(outlines)

	return None

if __name__ == '__main__':
	
# 	annotate_file_with_rep_timing_and_cpg('../../Mitra_etal_SFARI_SSC_denovo_TRs_Nov2020_reptiming_repeat_content.csv',
# 											'../allele_lengths.txt',
# 											'./test_rsdnm_annotation.txt')
	annotate_file_with_rep_timing_and_cpg(snakemake.input[0], snakemake.input[1], snakemake.output[0])
	
