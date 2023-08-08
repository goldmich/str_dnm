#!/usr/bin/env python

import pysam
import pdb
import re
import numpy as np
import os

locus_statistics = {}

infiles = snakemake.input[:22]
for f in infiles:
	with open(f, 'r') as open_infile:
		for l in open_infile.readlines():
			l_split = l.rstrip('\n').split('\t')
			if l_split[0] not in locus_statistics:
				locus_statistics[l_split[0]] = 0
			locus_statistics[l_split[0]] += int(l_split[1])

outfile = snakemake.output[0]
with open(outfile, 'w') as open_outfile:
	for ru in locus_statistics:
		open_outfile.write('%s\t%i\n' % (ru, locus_statistics[ru]))