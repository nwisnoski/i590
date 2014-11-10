#! /usr/bin/env python

usage = """
hw5 - Nathan Wisnoski

Takes a .fasta file and calculates LD between all pairs of SNPs.
Reports the number of SNPs, number of chromosomes, and value of r^2

"""

import sys
import itertools

def fasta_to_dict(in_file):
	# create empty list and dictionary to hold seqs
	in_seqs = {}
	current_seq = ""

	# cycle through lines in file and add to dictionary[ID] -> "sequence"
	for each in in_file:
		each = each.strip()
		if each.startswith(">"): # sequence ID
			current_seq = each
			in_seqs[current_seq] = "" # create dictionary key
		elif each == '///':
			pass
		else:
			in_seqs[current_seq] += each # add sequence to dictionary
	return in_seqs
	
def seqs_by_position(in_seqs):
	seqs = {}
	position = 0 # index through each position to setup the list
	
	# initialize dictionary of nucleotides at each locus
	for each in in_seqs.values():
		position = 0 # index from the beginning again
		for base in each:
			seqs[position] = []
			position += 1
	
	# add values to the dictionary
	for each in in_seqs.values():
		position = 0 # index from the beginning again
		for base in each:
			seqs[position] += base
			position += 1
	return seqs

def find_snps(in_seqs):
	nucleotides = set(['A','C','T','G'])
	snp_positions = []
	# identify SNP positions and calculate allele freqs, h, and fst
	for position in in_seqs.keys():
		set_position = set(in_seqs[position])
		if len(set_position & nucleotides) == 2:
			snp_positions.append(position)
	return snp_positions

def r_squared(seqs, snp_positions):
	nucleotides = ('A','C','T','G')
	r2 = []
	snp_pairs = itertools.combinations(snp_positions, 2)

	# cycle through all pairwise comparisons 
	for pair in snp_pairs:
		snp1_pos, snp2_pos = pair[0], pair[1]
		snp1, snp2 = seqs[snp1_pos], seqs[snp2_pos]
		# print snp1, snp2
		snp1_alleles = {}
		snp2_alleles = {}

		haplotypes = {}
		length1, length2 = 0, 0
		for each in snp1:
			if each != '-':
				length1 += 1
		for each in snp2:
			if each != '-':
				length2 += 1
		length = min(length1, length2)

		for position in range(length):
			haplotypes[position] = [snp1[position], snp2[position]]
		unique_haps = {}
		denominator = 0
		
		allele_freqs1 = {}
		allele_freqs2 = {}

		for hap in haplotypes.values():
			hap_str = ''.join(hap)
			if hap[0] and hap[1] in "ACTG":
				unique_haps[hap_str] = haplotypes.values().count(hap)/float(length)
				a1, a2 = hap[0], hap[1]
				if snp1.count(a1) > 0:
					allele_freqs1[a1] = snp1.count(a1)/float(length)
				if snp1.count(a2) > 0:
					allele_freqs1[a2] = snp1.count(a2)/float(length)
				if snp2.count(a1) > 0:
					allele_freqs2[a1] = snp2.count(a1)/float(length)
				if snp2.count(a2) > 0:
					allele_freqs2[a2] = snp2.count(a2)/float(length)
			elif hap[0] or hap[1] == '-':
				pass

		hap = unique_haps.keys()[1]
		a, b = hap[0], hap[1]
		ab = unique_haps[hap]
		a1 = allele_freqs1[a]
		b1 = 1 - a1
		a2 = allele_freqs2[b]
		b2 = 1 - a2
		D = ab - (a1 * a2)
		
		r2 = (D * D) / (a1 * b1 * a2 * b2)
		

		print str(snp1_pos + 1) + " and " + str(snp2_pos + 1) + ":\t r2=" + str(r2)
			
if len(sys.argv) < 2:
	print usage

else:
	# read and open filename supplied via command line
	in_file_name = sys.argv[1]
	print "Current file: " + in_file_name
	in_file = open(in_file_name, "r")
	
	seqs = fasta_to_dict(in_file)
	print "Number of Chromosomes: " + str(len(seqs.keys()))
	seqs = seqs_by_position(seqs)
	
	snp_positions = find_snps(seqs)
	print "Number of SNPs: " + str(len(snp_positions))

	print "LD:"
	r_squared(seqs, snp_positions)

	in_file.close()
	
	