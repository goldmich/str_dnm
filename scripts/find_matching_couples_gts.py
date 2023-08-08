#!/usr/bin/env python

import pysam
import pdb
import re
import numpy as np
import os
import csv
import random

dnm_file = '/net/harris/vol1/project/simons_simplex/Mitra_etal_SFARI_SSC_denovo_TRs_Nov2020_reptiming_repeat_content_asd_coding_regions.csv'
ped_file = '/net/harris/vol1/project/simons_simplex/ssc_copy.ped'
autosomes = ['chr%d' % (x) for x in range(1, 23)]

vcf_filename = snakemake.input[0]
# vcf_filename = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase3_1_10.sorted.filtered.vcf.gz'
curr_chr = 'chr' + re.search('(\d+)\.sorted', vcf_filename).group(1)

# Getting couples, saving family id this time
couples = []
with open(ped_file, 'r') as ped_in:
	fam_dict = {}
	for l in ped_in:
		l_split = l.rstrip('\n').split(' ')
		if l_split[1].endswith(('fa', 'mo')):
			if l_split[0] not in fam_dict:
				fam_dict[l_split[0]] = []
			fam_dict[l_split[0]].append(l_split[6])
	couples = [tuple([x[0]] + x[1]) for x in fam_dict.items()]

'''
read in dnms, sort into chromosome and position and save line info
for chr:
	open up relevant vcfs
	for dnm in chr:
		compile list of families (including this one) that share mutation here with same parental genotype (if parental genotypes are swapped that is okay)
		for each vcf:
			at the slice that includes the site of interest, are there any couples who match the gt requirements? if so, save their info
'''

dnm_dict = {}
with open(dnm_file) as open_dnm_file:
	open_csv = csv.reader(open_dnm_file, delimiter = ',', quotechar='"')
	next(open_csv)
	for l in open_csv:
		if l[0] != curr_chr:
			continue
		if l[1] not in dnm_dict:
			dnm_dict[l[1]] = {}
		# silly diy hashing bc sets aren't hashable in dict (can't be keys)
		parental_gts = '|'.join([min([l[12], l[13]]), max([l[12], l[13]])])
		if parental_gts not in dnm_dict[l[1]]:
			dnm_dict[l[1]][parental_gts] = [[], []]
		dnm_dict[l[1]][parental_gts][0].append(l[4])
		if l[8] not in dnm_dict[l[1]][parental_gts][1]:
			dnm_dict[l[1]][parental_gts][1].append(l[8])

vcf_stem = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp'
vcf_suffix = 'sorted.filtered.vcf.gz'

outlist = []
with pysam.VariantFile(vcf_filename) as vcf_in:
	couples_subset = [x for x in couples if all((x[1] in list(vcf_in.header.samples), x[2] in list(vcf_in.header.samples)))]
	for site in dnm_dict:
		for var in vcf_in.fetch(curr_chr, int(site), int(site) + 1):
# 			if random.random() < 0.01:
# 				pdb.set_trace()
			curr_period = var.info['PERIOD']
			allele_n_repeats = [int(len(x) / curr_period) for x in var.alleles]
			for gt in dnm_dict[site]:
				for family, s1, s2 in couples_subset:
					curr_gts = [var.samples[x]['GT'] for x in [s1, s2]]
					if any([len(x) < 2 for x in curr_gts]):
						continue
					gts_str = [','.join([str(allele_n_repeats[y]) for y in x]) for x in curr_gts]
					# yes i know that i'm min/maxing strings!! i do it above to save time and space
					curr_hash_gts = '|'.join([min([gts_str[0], gts_str[1]]), max([gts_str[0], gts_str[1]])])
					if curr_hash_gts == gt and family not in dnm_dict[site][gt][0]:
						max_gt = max([int(y) for x in curr_hash_gts.split('|') for y in x.split(',')] + 
									 [int(x) for x in dnm_dict[site][gt][1]])
						outlist.append([curr_chr, site, gt, ','.join(dnm_dict[site][gt][1]), 
										var.ref[:var.info['PERIOD']], str(max_gt), family, s1, s2])

# outfile = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase%s_%i.couples_for_str_distrib.txt' % (curr_phase, curr_chr)
outfile = snakemake.output[0]
with open(outfile, 'w')	as outfile_open:
# 	outfile_open.write('%s\t%s\t%s\t%s\t%s\n' % ('chr', 'pos', 'mut_family_gt', 'de_novo_alleles', 'repeat_unit', 'max_allele_length', 'family', 'id1', 'id2'))
	outfile_open.write('\n'.join(['\t'.join(x) for x in outlist]))
	outfile_open.write('\n')

