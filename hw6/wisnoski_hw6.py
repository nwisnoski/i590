#! /usr/bin/env python

usage = """
Homework 6 - Nathan Wisnoski
"""

import sys

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
def split_pop(in_dict):
	pop = []
	outgroup = []
	for each in in_dict.keys():
		if each.startswith(">WT"): # if this is from outgroup
			outgroup.append(in_dict[each])
		else: # if this is from pop A
			pop.append(in_dict[each])
	
	return pop, outgroup #returns separate lists of each population
def codons(in_seqs):
	current_seqs = []
	for each in in_seqs:
		# split sequences into groups of 3 - codons
		each = each.replace('T','U')
		current_seqs.append([each[i:i+3] for i in range(0, len(each), 3)])
	codons = {}
	for seq in current_seqs:
		position = 0
		for codon in seq:
			codons[position] = []
			position += 1
	for seq in current_seqs:
		position = 0
		for codon in seq:
			codons[position].append(codon)
			position += 1
	return codons
def differences(population, outgroup):
	bases = ['U', 'C', 'A', 'G']
	codons = [a+b+c for a in bases for b in bases for c in bases]
	amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
	codon_table = dict(zip(codons, amino_acids))

	Pn, Ps, Dn, Ds = 0,0,0,0

	for each in population.keys():
		codons = population[each]
		out = outgroup[each]
		if len(set(codons)) > 1: # there's a snp here
			poly = list(set(codons))
			if codon_table[poly[0]] == codon_table[poly[1]]: # syn poly
				Ps += 1
			else: # n.s. poly
				print codons
				Pn += 1
		if len(set(codons + out)) > 1:
			if out[0] not in codons:
				diffs = list(set(codons + out))
				if codon_table[diffs[0]] == codon_table[diffs[1]]: # syn diff
					Ds += 1
				else: # ns diff
					Dn += 1
	
	print 'Pn: \t' + str(Pn)
	print 'Ps: \t' + str(Ps)
	print 'Dn: \t' + str(Dn)
	print 'Ds: \t' + str(Ds)
	
if len(sys.argv) < 2:
	print usage
else:
	filename = sys.argv[1]
	in_file = open(filename, 'r')

	# parse fasta file into two lists by population
	seqs = fasta_to_dict(in_file)
	pop, outgroup = split_pop(seqs)
	pop = codons(pop)
	outgroup = codons(outgroup)
	differences(pop, outgroup)
	in_file.close()