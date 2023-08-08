#!/usr/bin/env python

import pysam
import pdb
import re
import numpy as np
import os

vcf_filename = snakemake.input[0]
# vcf_filename = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase3_1_10.sorted.filtered.vcf.gz'
chr = 'chr' + re.search('(\d+)\.sorted', vcf_filename).group(1)
samples = []
couple_dict = {}
gt_dict = {}
ped_file = '/net/harris/vol1/project/simons_simplex/ssc_copy.ped'

# Getting couples
couples = []
with open(ped_file, 'r') as ped_in:
	fam_dict = {}
	for l in ped_in:
		l_split = l.rstrip('\n').split(' ')
		if l_split[1].endswith(('fa', 'mo')):
			if l_split[0] not in fam_dict:
				fam_dict[l_split[0]] = []
			fam_dict[l_split[0]].append(l_split[6])
	couples = [tuple(x) for x in fam_dict.values()]

# Read in str dnm file here, save loci from chr for which we need allele info
dnm_filename = snakemake.input[1]
# dnm_filename = '/net/harris/vol1/project/simons_simplex/Mitra_etal_SFARI_SSC_denovo_TRs_Nov2020_reptiming_repeat_content.csv'
allele_length_dict = {}
with open(dnm_filename, 'r') as dnm_file_in:
	header = dnm_file_in.readline()
	for l in dnm_file_in.readlines():
		l_split = l.rstrip('\n').split(',')
		if l_split[0] != chr:
			continue
		# I don't care if loci get added more than once
		allele_length_dict[(l_split[0], int(l_split[1]))] = np.zeros((100), dtype = int)

# Check to make sure that either (a) there are fixed sites here or (b) segregating sites is an assumption we're okay with
# I need to test this a bit
# I can grab the repeat information from the vcf record!!!!!
with pysam.VariantFile(vcf_filename) as vcf_in:
	samples = list(vcf_in.header.samples)
	couples_subset = [(x1, x2) for x1, x2 in couples if all((x1 in samples, x2 in samples))]
	couple_dict = {x: np.zeros((5, 3), dtype = int) for x in couples_subset}
	counter = 0
	for var in vcf_in.fetch():
# 		counter += 1
# 		if counter > 1000:
# 			break
		curr_period = var.info['PERIOD']
		allele_n_repeats = [int(len(x) / curr_period) for x in var.alleles]
		# I'm defining the longest repeat unit as >=5; i can do this here because i already calculated n repeats above
		if curr_period >= 5:
			curr_period = 4
		# filters
		if not 'PASS' in var.filter:
			continue
		curr_locus = (var.chrom, var.pos)
		for s1, s2 in couples_subset:
			if len(var.samples[s1]['GT']) < 2 or len(var.samples[s2]['GT']) < 2:
				continue
			n_hets = sum([var.samples[s]['GT'][0] != var.samples[s]['GT'][1] for s in [s1, s2]])
			couple_dict[(s1, s2)][curr_period, n_hets] += 1
		if curr_locus in allele_length_dict:
			parental_alleles = [x for couple in couples_subset for s in couple for x in var.samples[s]['GT'] if len(var.samples[s]['GT']) > 1]
			parental_alleles_set = [allele_n_repeats[x] for x in sorted(list(set(parental_alleles)))]
			# i checked that the maximum # repeat units for an allele is 54 (49 in parents), so array out to 100 is fine
			try:
				allele_length_dict[curr_locus][parental_alleles_set] = 1
			except IndexError:
				allele_length_dict[curr_locus][[x for x in parental_alleles_set if x < 100]]

with open(snakemake.output[0], 'w') as open_couples_gt_outfile:
	header_list = [str(x) + "_" + str(y) for x in range(1, 5) for y in ['0_het', '1_het', '2_het']]
	open_couples_gt_outfile.write('sample1\tsample2\t' + '\t'.join(header_list) + '\n')
	for s1, s2 in couple_dict:
		open_couples_gt_outfile.write(s1 + '\t' + s2 + '\t')
		open_couples_gt_outfile.write('\t'.join([str(x) for x in np.reshape(couple_dict[(s1, s2)][1:,], -1)]) + '\n')

with open(snakemake.output[1], 'w') as open_allele_length_outfile:
# with open('test_allele_length_output.txt', 'w') as open_allele_length_outfile:
	open_allele_length_outfile.write('chrom\tpos\t' + '\t'.join([str(x) for x in range(100)]) + '\n')
	for chr, pos in sorted(allele_length_dict.keys(), key=lambda x: x[1]):
		open_allele_length_outfile.write(chr + '\t' + str(pos) + '\t' + '\t'.join(str(x) for x in allele_length_dict[(chr, pos)]) + '\n')
