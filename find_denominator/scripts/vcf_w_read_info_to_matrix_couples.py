#!/usr/bin/env python

import random
import pysam
import pdb
import re
import numpy as np
import os
import csv
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file")
ap.add_argument("-v", "--vcf", required=True, help="input file")
ap.add_argument("-f", "--family", required=True, help="family")
ap.add_argument("-c", "--chrom", required=True, help="chromosome")
ap.add_argument("-m", "--matrix_out", required=True, help="output file")
ap.add_argument("-d", "--dnms_out", required=True, help="output file")
args = ap.parse_args()

vcf_filename = args.vcf
# vcf_filename = '/net/eichler/vol26/projects/denovo_variation/nobackups/STRs/denominator/test/12220_chr22_uncallable.txt'
curr_family = args.family
curr_chr =  args.chrom
couple_matrix = np.zeros((5, 3, 2, 4), dtype = int)
parents = []

dnm_filename = args.input
# dnm_filename = '/net/eichler/vol26/projects/denovo_variation/nobackups/STRs/filtering/validated_denovo_STRs.csv'

family_dnm_lines = {}
with open(dnm_filename) as open_dnm_file:
	open_csv = csv.reader(open_dnm_file, delimiter = ',', quotechar='"')
	header = next(open_csv)
	curr_chr_prefix = 'chr' + curr_chr
	for l in open_csv:
		if l[0] != curr_chr_prefix or l[4] != curr_family:
			continue
		curr_pos = int(l[1])
		if curr_pos not in family_dnm_lines: # accounting for the hopefully very rare case where there are two DNMs reported in the kids at the same site, one in each
			family_dnm_lines[curr_pos] = []
		# just in case there's a dnm whose site gets filtered before adding discoverability info
		family_dnm_lines[curr_pos].append(l + [-1])


with pysam.VariantFile(vcf_filename) as vcf_in:
	parents = [x for x in vcf_in.header.samples]
	# repeat length, # hets, AT-only vs GC-containing, mutsize=1,2,3
	couple_matrix = np.zeros((5, 3, 2, 4), dtype = int)
	counter = 0
	for var in vcf_in.fetch():
# 		if random.random() < 0.001:
# 			pdb.set_trace()
		curr_period = var.info['PERIOD']
		allele_n_repeats = [int(len(x) / curr_period) for x in var.alleles]
		if curr_period >= 5:
			curr_period = 4
		# filters
		if not 'PASS' in var.filter:
			continue
		curr_repeat = var.ref[:var.info['PERIOD']]
		contains_gc = int('G' in curr_repeat or 'C' in curr_repeat)
		# determining heterozygosity
		if len(var.samples[parents[0]]['GT']) < 2 or len(var.samples[parents[1]]['GT']) < 2:
			continue
		n_hets = sum([var.samples[s]['GT'][0] != var.samples[s]['GT'][1] for s in parents])
		
		# possible mut sizes
		curr_seg_alleles = [allele_n_repeats[y] for s in parents for y in var.samples[s]['GT']]
		possible_mut_sizes = [0 for _ in range(4)]
		for mut_size in range(1, 4):
			possible_mut_sizes[mut_size] = int(all([(x - mut_size) not in var.info['UNCALLABLE_DNMS'] for x in curr_seg_alleles] + 
											   	   [(x + mut_size) not in var.info['UNCALLABLE_DNMS'] for x in curr_seg_alleles]))
		# fill in a matrix here
		couple_matrix[curr_period, n_hets, contains_gc] += possible_mut_sizes
		
		# add filter info to dnm file
		if var.pos in family_dnm_lines:
			for i in range(len(family_dnm_lines[var.pos])):
				if abs(int(family_dnm_lines[var.pos][i][9])) > 3:
					family_dnm_lines[var.pos][i][-1] = -1
				else:
					family_dnm_lines[var.pos][i][-1] = possible_mut_sizes[abs(int(family_dnm_lines[var.pos][i][9]))]

output_matrix_file = args.matrix_out
# output_matrix_file = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/%s.%s.parents_callable_dnms_matrix.txt' % (curr_family, curr_chr)

with open(output_matrix_file, 'w') as open_couples_gt_outfile:
	open_couples_gt_outfile.write('\t'.join(['sample1', 'sample2', 'ru_len', 'n_hets', 'contains_gc', 'dnm_size', 'count']) + '\n')
	for x in range(1, 5):
		for y in range(3):
			for z in range(2):
				for w in range(1, 4):
					open_couples_gt_outfile.write('\t'.join([str(x) for x in [parents[0], parents[1], x, y, z, w, couple_matrix[x, y, z, w]]]) + '\n')

output_dnm_lines_file = args.dnms_out
# output_dnm_lines_file = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/%s.%s.dnms_w_filter_info.csv' % (curr_family, curr_chr)
with open(output_dnm_lines_file, 'w') as open_dnm_lines_file:
	dnm_writer = csv.writer(open_dnm_lines_file, delimiter = ',', quotechar = '"')
	dnm_writer.writerow(header + ['all_dnms_discoverable'])
	for i in sorted(family_dnm_lines.keys()):
		dnm_writer.writerows(family_dnm_lines[i])
