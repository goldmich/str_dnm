#!/usr/bin/env python

import pysam
import pdb
import re
import numpy as np
import os
import csv
import gzip


# parents_reads_file = '/net/eichler/vol26/projects/denovo_variation/nobackups/new_quads/STR/distribution/new_sites/family_reads/13509.parents_reads_padding_2.txt'
parents_reads_file = snakemake.input[0]

site_dict = {}
with open(parents_reads_file, 'r') as open_file:
	open_file.readline() # header	
	for line in open_file.readlines():
		l_split = line.rstrip('\n').split()
		curr_hash = '_'.join(l_split[1:3])
		if curr_hash not in site_dict:
			site_dict[curr_hash] = l_split[4:7] + [0]
		site_dict[curr_hash][3] += int(l_split[8] in l_split[6].split(','))


# outfile = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/13509.parents_reads_distribution_padding_2.txt'
outfile = snakemake.output[0]
with open(outfile, 'w') as open_outfile:
	open_outfile.write('\n'.join(['\t'.join(site.split('_') + [str(x) for x in site_dict[site]]) for site in site_dict]))
	open_outfile.write('\n')