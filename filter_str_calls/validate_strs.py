import pysam
import argparse
from tqdm import tqdm
from collections import Counter
import pandas as pd

def get_bam(sample, manifest):
	"""
	Retrieves the path to the sample's bam file.
	:param sample: Sample to find bam file for.
	:returns: AlignmentFile object for the sample's bam.
	"""
	row = manifest[manifest['sample_id'] == sample]
	bam_path = row.bam.values[0]
	bam = pysam.AlignmentFile(bam_path, "rc")
	return bam

def get_family_bams(manifest, sample_list, family_id):
	"""
	Reads in all the alignment files (bams or crams) for a family.
	:param manifest: Pandas DataFrame with ssc_id, sample_id, and path to bam or cram file for each sample.
	:param sample_list: List of IDs for each family member.
	:param family_id: Integer of the unique ID for each family.
	:return: lLst of (name of member, AlignmentFile for bam) for mom, dad, pro, and sib.
	"""
	pro_id = [x for x in sample_list if 'p' in x][0]
	sib_id = [x for x in sample_list if 's' in x][0]

	mom = ("mom", get_bam("%s.mo" % family_id, manifest))
	dad = ("dad", get_bam("%s.fa" % family_id, manifest))
	pro = (pro_id, get_bam(pro_id, manifest))
	sib = (sib_id, get_bam(sib_id, manifest))

	return [mom, dad, pro, sib]


def get_longest_genotype(variant_data):
	"""
	Find the length of longest genotype.
	:param variant_data: DataFrame row for a variant, with child_gt, mat_gt, and pat_gt columns.
	:return: Integer of the length of longest genotype (STR motif count multiplied by motif length).
	"""
	motif_length = int(variant_data['len_repeat_unit'])
	genotypes = []
	for gt in ["child_gt", "mat_gt", "pat_gt"]:
		genotypes += [int(x) for x in variant_data[gt].split(',')]
	return max(genotypes) * motif_length

def get_denovo_genotype(variant_data):
	"""
	Get the number of motif repeats in the de novo genotype.
	:param variant_data: DataFrame row for a variant, with child_gt, mat_gt, and pat_gt columns.
	:return: List of integer number of repeat units in de novo genotype, and parent of origin if known.
	"""
	child_genotypes = set(variant_data['child_gt'].split(','))
	parental_genotypes = set(variant_data['pat_gt'].split(',') + variant_data['mat_gt'].split(','))
	denovo_genotype = child_genotypes - parental_genotypes

	if len(denovo_genotype) != 0:
		return [int(x) for x in denovo_genotype], None

	# If there are no alleles in the child that are missing in the parents, that means that both of the child's
	# Alleles have come from one parent. Since we don't know which one is de novo in this case, return both alleles
	# And specify which parent the de novo allele is from.
	if child_genotypes - set(variant_data['pat_gt'].split(',')) :
		return [int(x) for x in child_genotypes], 'dad'
	elif child_genotypes - set(variant_data['mat_gt'].split(',')) :
		return [int(x) for x in child_genotypes], 'mom'
	else:
		return [int(x) for x in child_genotypes], None

def get_soft_clipping(read):
	"""
	Gets the number of soft clipped bases on either end of the read.
	:param read: AlignedSegment of sequencing read.
	:return: Integer lengths of leading and tailing soft clipped bases.
	"""
	cigar = read.cigartuples

	if cigar[0][0] == 4:
		leading_clip = cigar[0][1]
	else:
		leading_clip = 0

	if cigar[-1][0] == 4:
		tailing_clip = cigar[-1][1]
	else:
		tailing_clip = 0

	return leading_clip, tailing_clip

def check_cigar_start(read, str_start, longest_gt, padding):
	"""
	Checks to see if their are bases that do not match the reference before the STR start.
	:param read: AlignedSegment of sequencing read.
	:param str_start: Integer naive prediction of STR starting position in read.
	:param longest_gt: Integer length of longest observed STR genotype.
	:param padding: Integer length of padding to add to either side of STR sequence.
	:return: Integer STR start and end positions adjusted to accommodate for inserted/deleted bases earlier in the read.
	"""
	start_pos = 0
	adjusted_str_start = str_start

	for operation, length in read.cigartuples:
		# Once we reach the start, we no longer need to adjust its location
		if adjusted_str_start - 1 < start_pos + length :
			left_boundary = adjusted_str_start - padding
			right_boundary = adjusted_str_start + longest_gt + padding
			# If there's an insertion or deletion that spans the starting coordinate, include those bases in the sequence
			if operation in (1,2,3):
				left_boundary -= length
				right_boundary += length
			break
		# If there is an insertion (1) or soft clip (4), the position of the STR will be later in the read
		if operation in (1, 4):
			adjusted_str_start += length
		# If there is a deletion (2) or reference skip (3), the position of the STR will be earlier in the read
		elif operation in (2, 3):
			adjusted_str_start -= length

		start_pos += length

	return left_boundary, right_boundary

def get_read_sequence(read, str_coordinate, longest_gt, padding):
	"""
	Checks that reads pass mapq and span the STR, and pulls out the STR sequence with padding.
	:param read: AlignedSegment of sequencing read.
	:param str_coordinate: Integer for STR starting coordinate in reference.
	:param longest_gt: Number of bases in longest observed STR genotype.
	:param padding: Number of bases to add on either side of the STR.
	:return: String for STR sequence, with length of longest_gt + 2*padding, or None if read does not pass filters.
	"""
	# Filter on some of the same basic features that GangSTR uses
	if read.mapping_quality != 60 or read.is_secondary or read.is_supplementary:
		return

	# Find the start and end coordinates of the STR sequence with padding, based on the read alignment and CIGAR string
	naive_start = str_coordinate - read.reference_start - 1
	left, right = check_cigar_start(read, naive_start, longest_gt, padding)
	left_clip, right_clip = get_soft_clipping(read)

	# Check that the read spans the STR sequence and padding, excluding soft-clipped bases
	if left < left_clip or right > read.query_length - right_clip:
		return

	# print(read.query_name)
	# print(read.query_sequence[left:right])
	return read.query_sequence[left:right]

def check_interruption(sequence, start_pos, padding, prev_sequence, next_sequence, str_length):
	"""
	Check to see if the sequence following the STR matches the expected reference sequence.
	:param sequence: Read sequence.
	:param start_pos: Position of the base that interrupts the STR. Should be STR end + 1 for uninterrupted STR.
	:param padding: Padding around motif sequence.
	:param prev_sequence: The expected reference sequence preceding the STR.
	:param next_sequence: The expected reference sequence following the STR.
	:return: True if the following sequence matches expectation, False if it does not.
	"""
	start, end = True, True

	start_sequence = sequence[start_pos-padding-str_length + 1: start_pos-str_length]
	end_sequence = sequence[start_pos: start_pos + padding - 1]

	if len(start_sequence) != padding - 1:
		start = False
	if len(end_sequence) != padding - 1:
		end = False

	if start or end:
		start_differences = sum(c1 != c2 for c1, c2 in zip(start_sequence, prev_sequence))
		end_differences = sum(c1 != c2 for c1, c2 in zip(end_sequence, next_sequence))

		# Allow for up to 1 difference between expected and observed to account for base calling error around STR.
		# If there is an error in the first motif copy, the first motif will be included in the preceding sequence, so you should never get a match.
		if start_differences > 1:
			start = False
		# If it is an interrupted repeat, the base with the error will always be the first base, so you should never get a match with the subsequent sequence, regardless of similarity.
		if end_differences > 1:
			end = False

	return start, end

def find_motif_count(sequence, motif, longest_gt, padding, lead_seq, tail_seq):
	"""
	Counts the number of contiguous repeat motifs in a read.
	:param sequence: String of the sequence from a read.
	:param motif: String of the STR repeat motif.
	:param padding: Integer distance from the STR reference start coordinate to the start of the sequence.
	:param lead_seq: The 10 bp of sequence prior to the STR.
	:param tail_seq: The 10 bp of sequence after the STR.
	:return: Integer count of contiguous repeats of STR motif.
	"""
	start = end = 0
	count = 0

	# extra_padding = len(sequence) - (longest_gt + 2 * padding)
	extra_padding = 1

	while start <= len(sequence) - len(motif) :
		end += 1
		if sequence[end - 1] == motif[end - start - 1]:
			if end - start == len(motif):
				# Found an instance of the motif, count it.
				count += 1
				start = end
		else:
			# This is intended to address any alignment issues, so if an STR begins at a coordinate before the reference, those motifs will still be counted.
			if end - 1 <= padding + extra_padding:
				# Motif in the sequence before STR coordinates begin, and it is not contiguous with the STR.
				start = end = start + 1
				count = 0
			elif count == 0:
				# Keep going because you have not found the STR start yet.
				start = end = start + 1
			elif count != 0:
				# If there is a base error in the STR sequence, we don't want to count the length of that STR.
				start_ok, end_ok = check_interruption(sequence, start, padding, lead_seq, tail_seq, count*len(motif))
				if start_ok and end_ok:
					return count
				break
	return -1

def check_reads(bam_info, chrom, pos, longest_gt, motif, lead_seq, tail_seq, flank_length):
	"""
	Iterate through all reads aligned to STR coordinate and find the number of repeated motifs in each read.
	:param bam_info: Tuple of sample name and AlignmentFile object.
	:param chrom: Chromosome name.
	:param pos: STR start coordinate.
	:param longest_gt: Length of the longest observed STR genotype.
	:param motif: Repeated sequence in STR.
	:param lead_seq: The of sequence prior to the STR.
	:param tail_seq: The of sequence after the STR.
	:param flank_length: Length of the lead and tail seqs.
	:return: List of number of motifs observed in each read.
	"""
	sample_name, bam = bam_info
	motif_counts = []

	padding = flank_length + 1

	try:
		# print(sample_name)
		for aligned_segment in bam.fetch(chrom, pos - 1, pos):
			# print(aligned_segment.query_name)
			str_sequence = get_read_sequence(aligned_segment, pos, longest_gt, padding)
			# print(str_sequence)
			if str_sequence is not None:
				motif_ct = (find_motif_count(str_sequence, motif, longest_gt, padding, lead_seq, tail_seq))
				# print(motif_ct)
				if motif_ct != -1 :
					motif_counts.append(motif_ct)

		return motif_counts

	# if the CRAM is corrupted, return NA
	except OSError:
		return ['NA']

def get_denovo_motif_counts(family_reads, genotype, child):
	other_child = list(read_dict.keys() - ['mom', 'dad', denovo_child])[0]
	mom_motifs = family_reads['mom'].count(genotype)
	dad_motifs = family_reads['dad'].count(genotype)
	denovo_child_motifs = family_reads[child].count(genotype)
	other_child_motifs = family_reads[other_child].count(genotype)

	return mom_motifs, dad_motifs, denovo_child_motifs, other_child_motifs

def validate(mom_motif, dad_motif, child_motif, other_child_motif, threshold_allele_count):
	"""
	Decide whether the de novo mutation is true, inherited, or false based on the number of reads with the mutation in the parents and child.
	:param mom_motif: The number of reads with the de novo motif in the mom.
	:param dad_motif: The number of reads with the de novo motif in the dad.
	:param child_motif: The number of reads with the de novo motif in the child.
	:param threshold_allele_count: Threshold number of reads with the de novo allele to say that a parent definitely has it.
	:return: Validation status of de novo allele based on read data, and parent_of_origin if determined.
	"""
	if child_motif == 0:
		return 'false_positive', 'unknown'

	# The mom has enough reads with the de novo allele to confidently call it inherited
	if mom_motif >= threshold_allele_count :
		if dad_motif < threshold_allele_count :
			return 'inherited_from_mom', 'mom'
		else:
			return 'inherited_from_mom_or_dad', 'mom_or_dad'

	# The dad has enough reads with the de novo allele to confidently call it inherited
	if dad_motif >= threshold_allele_count :
		if mom_motif < threshold_allele_count :
			return 'inherited_from_dad', 'dad'
		else:
			return 'inherited_from_mom_or_dad', 'mom_or_dad'

	# The mom has some reads with the de novo allele, but not enough to be confident - so we'll equivocate!
	if 0 < mom_motif < threshold_allele_count :
		if dad_motif == 0 :
			return 'potentially_inherited_from_mom', 'mom'
		else:
			return 'potentially_inherited_from_mom_or_dad', 'mom_or_dad'

	# The dad has some reads with the de novo allele, but not enough to be confident - so we'll equivocate!
	if 0 < dad_motif < threshold_allele_count:
		if mom_motif == 0 :
			return 'potentially_inherited_from_dad', 'dad'
		else:
			return 'potentially_inherited_from_mom_or_dad', 'mom_or_dad'

	if other_child_motif > 0:
		return 'inherited', 'unknown'

	# The parents have no reads with the de novo allele, so it's as real as it gets
	if mom_motif == 0 and dad_motif == 0:
		return 'true_de_novo', 'unknown'

def filter_motif_counts(family_reads, denovo_genotype, naive_parent, denovo_child, threshold):
	"""
	Check that the de novo motif is present in the child and not present in either parent.
	:param family_reads: Dictionary with keys of sample name and values of list of number of STR motifs in each read.
	:param denovo_genotype: Integer of number of motifs in the de novo genotype.
	:param naive_parent: Naive prediction of parent of origin, None if not determined.
	:param denovo_child: Sample name of child with the de novo mutation.
	:param threshold: Number of reads with de novo allele to establish inheritance.
	:return: Validation status of de novo mutation.
	"""
	# We know the de novo genotype, so it's a really simply validation
	if len(denovo_genotype) == 1:
		mom, dad, denovo_kid, other_kid = get_denovo_motif_counts(family_reads, denovo_genotype[0], denovo_child)
		validation, parent_of_origin = validate(mom, dad, denovo_kid, other_kid, threshold)
		return validation

	# The de novo genotype might be one of two alleles, and we can predict the parent of origin based on allele presence/absence
	if len(denovo_genotype) == 2 and naive_parent:
		mom_counts = []
		dad_counts = []
		# Check each child genotype to see if it is de novo or inherited
		for gt in denovo_genotype:
			mom, dad, denovo_kid, other_kid = get_denovo_motif_counts(family_reads, gt, denovo_child)
			# if at least one allele is observed in both parents, the other is likely inherited
			if mom != 0 and dad != 0:
				return 'inherited_from_mom_or_dad'
			else:
				mom_counts.append(mom)
				dad_counts.append(dad)

		if mom_counts[0] != 0:
			if dad_counts[1] != 0:
				return "inherited"
		if dad_counts[0] != 0:
			if mom_counts[1] != 0:
				return "inherited"

		return 'true_de_novo'

	return 'unknown'

def get_allele_count(motif_counts):
	"""
	Summarize motif counts to get the number of reads that support each observed allele.
	:param motif_counts: List of a motif repeat count for each spanning read.
	:return:  String of each allele and number of supporting reads separated by commas.
	"""
	allele_counts = Counter(motif_counts)
	allele_list = []
	for allele in allele_counts:
		allele_list.append(str(allele) + '=' + str(allele_counts[allele]))

	return ','.join(allele_list)

def get_sample_name(sample):
	"""
	Set sample name to pro or sib if sample ID is in SSC format, in order to make generic columns.
	:param sample: Sample ID.
	:return: Simplified sample name.
	"""
	if 'p' in sample:
		return 'pro'
	elif 's' in sample:
		return 'sib'
	else:
		return sample

def get_flanking_sequences(str_chrom, str_pos, motif_len, ref_len, reference, padding):
	"""
	Get the reference base pairs before/after the STR.
	:param str_chrom: Chromosome of the STR.
	:param str_pos: Starting coordinate of the STR.
	:param motif_len: Length of the STR motif.
	:param ref_len: Number of STR copies in the reference genome.
	:param reference: PySam object of reference genome.
	:param padding: Length of flanking sequence.
	:return: Strings of the 10 base pairs before and 10 base pairs after the STR.
	"""
	# Subtract 1 from the positions to correctly index bam with pysam
	start = str_pos - 1
	end = str_pos + (motif_len * ref_len) - 1
	prev_seq = reference.fetch(str_chrom, start - padding, start)
	next_seq = reference.fetch(str_chrom, end, end + padding)
	return next_seq.upper(), prev_seq.upper()


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
	ap.add_argument("-b", "--bams", required=True, help="Manifest of bam files")
	ap.add_argument("-t", "--threshold", required=False, help="Number of reads with de novo allele to establish inheritance",
					default=1, type=int)
	ap.add_argument("-o", "--output", required=False, help="Manifest of bam files", default="STR_validations.csv")
	ap.add_argument("-r", "--reference", required=False, help="Reference genome fasta")
	ap.add_argument("-f", "--flanking_length", required=False, help="Length of flanking sequence",
					default=5, type=int)
	args = ap.parse_args()

	bams = pd.read_csv(args.bams, sep='\t')
	variants = pd.read_csv(args.input)
	ref = pysam.FastaFile(args.reference)

	families = list(set(variants['family']))

	for family in tqdm(families, total=len(families)):
		family_variants = variants.loc[variants['family'] == family]
		samples = bams.loc[bams['sample_id'].str.contains(str(family))]['sample_id'].values

		# There are some families that I cannot find the crams for, so just move along
		if samples.size == 0:
			print("Missing manifest data for family", str(family))
			continue

		family_bams = get_family_bams(bams, samples, family)

		for idx, row in family_variants.iterrows():
			# print(row['chrom'], row['pos'])
			denovo_child = bams.loc[bams['ssc_id'] == row['child']]['sample_id'].values[0]
			ref_length = int(row['ref_allele_len'])
			longest_gt = get_longest_genotype(row)
			denovo_gt, parent_of_origin = get_denovo_genotype(row)
			next_seq, prev_seq = get_flanking_sequences(row['chrom'], row['pos'], row['len_repeat_unit'], ref_length, ref, args.flanking_length)

			if not denovo_gt:
				validation = "no_denovo_gt_listed"
				continue

			read_dict = {}

			for bam in family_bams:
				# print(bam[0])
				read_dict[bam[0]] = check_reads(bam, row['chrom'], row['pos'], longest_gt, row['repeat_unit'], prev_seq, next_seq, args.flanking_length)
				# Add the read counts for each allele to a column in the DataFrame.
				variants.loc[idx, "%s_allele_counts"%get_sample_name(bam[0])] = get_allele_count(read_dict[bam[0]])

			validation = filter_motif_counts(read_dict, denovo_gt, parent_of_origin, denovo_child, args.threshold)
			variants.loc[idx, 'validation'] = validation

	variants.to_csv(args.output, index = False)
