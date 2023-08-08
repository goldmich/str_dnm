#!/usr/bin/env python

import pdb
import re
import numpy as np
import csv
import os

def annotate_panel():
	
	panel_file = '/net/harris/vol1/project/simons_simplex/hg38_ver16.bed'
	reptiming_file = '/net/harris/vol1/project/simons_simplex/koren_reptiming_data/ESC.spar1e-16.PC24.hg38.consensus.bed'
	dnm_file = '/net/harris/vol1/project/simons_simplex/full_validation_padding_2.csv'
	autosomes = ['chr%i' % x for x in range(1, 23)]
	
	dnm_dict = {}
	
	with open(dnm_file) as open_dnm_file:
		dnm_reader = csv.reader(open_dnm_file, delimiter = ",")
		dnm_reader.__next__()
		for line in dnm_reader:
			curr_chr = line[0]
			# autosomes only, no homopolymers, must pass filters
			if curr_chr not in autosomes or line[30] != 'true_de_novo' or line[2] == '1':
				continue
			if curr_chr not in dnm_dict:
				dnm_dict[curr_chr] = {}
			line[1] = int(line[1])
			if line[1] not in dnm_dict[curr_chr]:
				dnm_dict[curr_chr][line[1]] = 0
			dnm_dict[curr_chr][line[1]] += 1
	
	# normally i would just loop through the reptiming and str files but they are unsorted
	reptiming_dict = {}
	with open(reptiming_file, 'r') as open_reptiming_file:
		for line in open_reptiming_file.readlines():
			line_split = line.rstrip('\n').split('\t')
			if line_split[0] not in reptiming_dict:
				reptiming_dict[line_split[0]] = []
			reptiming_dict[line_split[0]].append([int(line_split[1]), int(line_split[2]), line_split[3]])
	
	for chr in reptiming_dict:
		reptiming_dict[chr] = sorted(reptiming_dict[chr], key=lambda x: x[0])
	
	outlines = []
	reptiming_pointer = 0
	with open(panel_file, 'r') as open_panel_file:
		for line in open_panel_file:
			line_split = line.rstrip('\n').split('\t')
			# new chr?
			if line_split[0] != curr_chr:
				reptiming_pointer = 0
				curr_chr = line_split[0]
			# autosomes only
			if curr_chr not in autosomes:
				line_split += ['NA', 'NA']
				outlines.append(line_split)
				continue
			curr_pos = int(line_split[1])
			# advancing reptiming pointer to the first segment where the end is not before the curr_pos
			while reptiming_pointer < len(reptiming_dict[curr_chr]) and curr_pos > reptiming_dict[curr_chr][reptiming_pointer][1]:
				reptiming_pointer += 1
			# if curr_pos is in reptiming window, append
			if reptiming_pointer < len(reptiming_dict[curr_chr]) and curr_pos >= reptiming_dict[curr_chr][reptiming_pointer][0]:
				line_split.append(reptiming_dict[curr_chr][reptiming_pointer][2])
			else:
				line_split.append('NA')
			
			# Any dnms?
			if curr_pos in dnm_dict[line_split[0]]:
				line_split.append(dnm_dict[line_split[0]][curr_pos])
			else:
				line_split.append(0)
			outlines.append(line_split)
	
	outfile = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/hg38_ver16.reptiming.dnm_rate.bed'
	with open(outfile, 'w') as open_outfile:
		open_outfile.write('\n'.join(['\t'.join([str(x) for x in y]) for y in outlines]))
		open_outfile.write('\n')
				
	

if __name__ == '__main__':
	annotate_panel()