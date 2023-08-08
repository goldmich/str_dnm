#!/usr/bin/env python

import pysam
import pdb
import re
import numpy as np
import os

vcf_filename = snakemake.input[0]
# vcf_filename = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase3_1_10.sorted.filtered.vcf.gz'
vcf_stem = '/net/harris/vol1/project/simons_simplex/indiv_str_distribution/tmp/phase3_1_'
vcf_suffix = '.sorted.filtered.vcf.gz'
# autosomes = [str(x) for x in range(1, 23)]
locus_statistics = {}

with pysam.VariantFile(vcf_filename) as vcf_in:
	for var in vcf_in.fetch():
		if not 'PASS' in var.filter:
			continue
		curr_ru = var.ref[:var.info['PERIOD']]
		if curr_ru not in locus_statistics:
			locus_statistics[curr_ru] = 0
		locus_statistics[curr_ru] += 1

outfile = snakemake.output[0]
with open(outfile, 'w') as open_outfile:
	for ru in locus_statistics:
		open_outfile.write('%s\t%i\n' % (ru, locus_statistics[ru]))