#!/usr/bin/env python

import pysam
import pdb
import re
import numpy as np
import os

# vcf_filename = snakemake.input[0]
vcf_filename = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase3_1_10.sorted.filtered.vcf.gz'
chr = 'chr' + re.search('(\d+)\.sorted', vcf_filename).group(1)
samples = []
sample_dict = {}
gt_dict = {}
allele_length_dict = {}

# Check to make sure that either (a) there are fixed sites here or (b) segregating sites is an assumption we're okay with
# I need to test this a bit
# I can grab the repeat information from the vcf record!!!!!
with pysam.VariantFile(vcf_filename) as vcf_in:
	samples = list(vcf_in.header.samples)
	sample_dict = {x: np.zeros((5, 100), dtype = int) for x in samples}
	gt_dict = {x: np.zeros((5, 2), dtype = int) for x in samples}
	counter = 0
	for var in vcf_in.fetch():
# 		counter += 1
# 		if counter > 1000:
# 			break
		if not 'PASS' in var.filter:
			continue
		curr_period = var.info['PERIOD']
		allele_n_repeats = [int(len(x) / curr_period) for x in var.alleles]
		# I'm defining the longest repeat unit as >=5; i can do this here because i already calculated n repeats above
		if curr_period >= 5:
			curr_period = 4
		for s in samples:
			if len(var.samples[s]['GT']) < 2:
				continue
			# homozygous
			if var.samples[s]['GT'][0] == var.samples[s]['GT'][1]:
				gt_dict[s][curr_period, 0] += 1
			# heterozygous
			else:
				gt_dict[s][curr_period, 1] += 1
			for allele in var.samples[s]['GT']:
				# the first time i see an allele of length >= 100 i will set it to 100
				if allele_n_repeats[allele] > 99:
					allele_n_repeats[allele] = 99
				sample_dict[s][curr_period, allele_n_repeats[allele]] += 1

with open(snakemake.output[0], 'w') as open_outfile:
	header_list = [str(x) + "_" + str(y) for x in range(1, 5) for y in range(100)]
	open_outfile.write('sample\t' + '\t'.join(header_list) + '\n')
	for s in sample_dict:
		open_outfile.write(s + '\t')
		open_outfile.write('\t'.join([str(x) for x in np.reshape(sample_dict[s][1:,], -1)]) + '\n')

with open(snakemake.output[1], 'w') as open_gt_outfile:
	header_list = [str(x) + "_" + str(y) for x in range(1, 5) for y in ['homo', 'het']]
	open_gt_outfile.write('sample\t' + '\t'.join(header_list) + '\n')
	for s in gt_dict:
		open_gt_outfile.write(s + '\t')
		open_gt_outfile.write('\t'.join([str(x) for x in np.reshape(gt_dict[s][1:,], -1)]) + '\n')
