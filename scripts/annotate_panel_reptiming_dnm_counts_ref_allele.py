#!/usr/bin/env python

import pysam
import pdb
import re
import numpy as np
import os
import csv

vcf_filename = snakemake.input[0]
# vcf_filename = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase3_1_10.sorted.filtered.vcf.gz'
vcf_stem = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase3_1_'
vcf_suffix = '.sorted.filtered.vcf.gz'
curr_chr = 'chr' + os.path.basename(vcf_filename).split('.')[0].split('_')[-1]
# autosomes = [str(x) for x in range(1, 23)]

dnm_file = '/net/harris/vol1/project/simons_simplex/full_validation_padding_2.csv'

dnm_dict = {}
# Need to change this to be for a single chromosome only
with open(dnm_file) as open_dnm_file:
	dnm_reader = csv.reader(open_dnm_file, delimiter = ",")
	dnm_reader.__next__()
	for line in dnm_reader:
		# autosomes only, no homopolymers, must pass filters
		if line[0] != curr_chr  or line[30] != 'true_de_novo' or line[2] == '1':
			continue
		line[1] = int(line[1])
		if line[1] not in dnm_dict:
			dnm_dict[line[1]] = 0
		dnm_dict[line[1]] += 1

# normally i would just loop through the reptiming and str files but they are unsorted
reptiming_list = []
reptiming_file = '/net/harris/vol1/project/simons_simplex/koren_reptiming_data/ESC.spar1e-16.PC24.hg38.consensus.bed'
with open(reptiming_file, 'r') as open_reptiming_file:
	for line in open_reptiming_file.readlines():
		line_split = line.rstrip('\n').split('\t')
		if line_split[0] != curr_chr:
			continue
		reptiming_list.append([int(line_split[1]), int(line_split[2]), line_split[3]])

reptiming_list = sorted(reptiming_list, key=lambda x: x[0])

outlines = []
with pysam.VariantFile(vcf_filename) as vcf_in:
	reptiming_pointer = 0
	for var in vcf_in.fetch():
		# no homopolymers, variants that don't pass filters
		if not 'PASS' in var.filter or var.info['PERIOD'] == 1:
			continue
		curr_pos = var.pos
		curr_ru = var.ref[:var.info['PERIOD']]
		curr_ref_len = int(len(var.ref) / var.info['PERIOD'])
		curr_line = [var.chrom, var.pos, curr_ru, curr_ref_len]
		# advancing reptiming pointer to the first segment where the end is not before the curr_pos
		while reptiming_pointer < len(reptiming_list) and curr_pos > reptiming_list[reptiming_pointer][1]:
			reptiming_pointer += 1
		# if curr_pos is in reptiming window, append
		if reptiming_pointer < len(reptiming_list) and curr_pos >= reptiming_list[reptiming_pointer][0]:
			curr_line.append(reptiming_list[reptiming_pointer][2])
		else:
			curr_line.append('NA')
		
		# Any dnms?
		if curr_pos in dnm_dict:
			curr_line.append(dnm_dict[curr_pos])
		else:
			curr_line.append(0)
		outlines.append(curr_line)

outfile = snakemake.output[0]
# outfile = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/sites_annotated_reptiming_dnm_counts_ru.10.txt'
with open(outfile, 'w') as open_outfile:
	open_outfile.write('\n'.join(['\t'.join([str(x) for x in y]) for y in outlines]) + '\n')
