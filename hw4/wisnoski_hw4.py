#! /usr/bin/env python

usage = """
hw4 - Nathan Wisnoski

Takes a FASTA-formated file for two populations.
Calculates Fst between the two populations.

Usage: "wisnoski_hw2.py sequence.fa"
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

def split_pops(in_dict):
	pop_A = []
	pop_B = []
	for each in in_dict.keys():
		if each.startswith(">B"): # if this is from pop B
			pop_B.append(in_dict[each])
		else: # if this is from pop A
			pop_A.append(in_dict[each])
	
	return pop_A, pop_B #returns separate lists of each population

def seqs_by_position(population):
	seqs = {}
	position = 0 # index through each position to setup the list
	
	# initialize dictionary of nucleotides at each locus
	for sequence in population:
		position = 0 # index from the beginning again
		for base in sequence:
			seqs[position] = []
			position += 1
	
	# add values to the dictionary
	for sequence in population:
		position = 0 # index from the beginning again
		for base in sequence:
			seqs[position] += base
			position += 1
	return seqs #returns the sequences arranged by position

def fst(pop_A, pop_B):
	nucleotides = set(['A','C','T','G'])
	n_A = len(pop_A.values())
	n_B = len(pop_B.values())
	Fst = {}
	
	# combine pops to calculate total h
	pop_pooled = {}
	for key in pop_A.keys():
		pop_pooled[key] = pop_A[key] + pop_B[key]
	n_pooled = n_pop_A + n_pop_B
	
	# identify SNP positions and calculate allele freqs, h, and fst
	for position in pop_pooled.keys():
		set_position = set(pop_pooled[position])
		if len(set_position) > 1:
			if set_position <= nucleotides: # SNP here
				pA_A = pop_A[position].count('A') / float(len(pop_A[position]))
				pA_C = pop_A[position].count('C') / float(len(pop_A[position]))
				pA_T = pop_A[position].count('T') / float(len(pop_A[position]))
				pA_G = pop_A[position].count('G') / float(len(pop_A[position]))
				pB_A = pop_B[position].count('A') / float(len(pop_B[position]))
				pB_C = pop_B[position].count('C') / float(len(pop_B[position]))
				pB_T = pop_B[position].count('T') / float(len(pop_B[position]))
				pB_G = pop_B[position].count('G') / float(len(pop_B[position]))
				
				# lists containing allele frequencies at this position for pops A & B
				freq_A = [pA_A,pA_C,pA_T,pA_G]
				freq_B = [pB_A,pB_C,pB_T,pB_G]
				
				# square the terms for Pi^2 values in h calculations
				pA_squared = []
				for each in freq_A:
					pA_squared.append(float(each) * float(each))
				pB_squared = []
				for each in freq_B:
					pB_squared.append(float(each) * float(each))	
				pbar = []
				for each in range(4):
					pbar_i = (freq_A[each]*n_pop_A + freq_B[each]*n_pop_B)/n_pooled
					pbar.append(pbar_i * pbar_i)
				
				# calculate ht and hw at this position
				ht = n_pooled/(n_pooled-1) * (1-(sum(pbar)))
				hw = ((n_pop_A/(n_pop_A-1)*(1-sum(pA_squared))) + (n_pop_B/(n_pop_B-1)*(1-sum(pB_squared))))/2
				
				# calculate Fst and add to dictionary indexed by position
				fst = (ht - hw)/ht
				Fst[position] = fst
			# repeat the same process as above, but for positions with indels
			elif len(set_position & nucleotides) > 1: # SNP here
				indelsA = pop_A[position].count('-')
				indelsB = pop_B[position].count('-')
				pA_A = pop_A[position].count('A') / float(len(pop_A[position])-indelsA)
				pA_C = pop_A[position].count('C') / float(len(pop_A[position])-indelsA)
				pA_T = pop_A[position].count('T') / float(len(pop_A[position])-indelsA)
				pA_G = pop_A[position].count('G') / float(len(pop_A[position])-indelsA)
				pB_A = pop_B[position].count('A') / float(len(pop_B[position])-indelsB)
				pB_C = pop_B[position].count('C') / float(len(pop_B[position])-indelsB)
				pB_T = pop_B[position].count('T') / float(len(pop_B[position])-indelsB)
				pB_G = pop_B[position].count('G') / float(len(pop_B[position])-indelsB)
				freq_A = [pA_A,pA_C,pA_T,pA_G]
				freq_B = [pB_A,pB_C,pB_T,pB_G]
				
				pA_squared = []
				for each in freq_A:
					pA_squared.append(float(each) * float(each))
				pB_squared = []
				for each in freq_B:
					pB_squared.append(float(each) * float(each))	
				pbar = []
				for each in range(4):
					pbar_i = (freq_A[each] + freq_B[each])/2
					pbar.append(pbar_i * pbar_i)
				
				ht = n_pooled/(n_pooled-1) * (1-(sum(pbar)))
				hw = ((n_pop_A/(n_pop_A-1)*(1-sum(pA_squared))) + (n_pop_B/(n_pop_B-1)*(1-sum(pB_squared))))/2
				fst = (ht - hw)/ht
				Fst[position] = fst
	# print the output to the terminal
	for each in Fst.keys():
		if Fst[each] < 0:
			Fst[each] = 0.0
		print "Position " + str(each) + '\t' + 'Fst=' + str(Fst[each])
		
# checks number of arguments
if len(sys.argv) < 2:
	print usage
else:
	filename = sys.argv[1]
	in_file = open(filename, 'r')
	
	# parse fasta file into two lists by population
	seqs = fasta_to_dict(in_file)
	pop_A, pop_B = split_pops(seqs)
	
	# calculate population size
	n_pop_A, n_pop_B = len(pop_A), len(pop_B)
	
	# organize nucleotides by position
	pop_A = seqs_by_position(pop_A)
	pop_B = seqs_by_position(pop_B)
	
	# calculate Fst
	Fst = fst(pop_A, pop_B)
	
	in_file.close()
