#! /usr/bin/env python

usage = """
hw1 - Nathan Wisnoski

Takes a FASTA-formated file and reports:
- Number of SNPs in the sample
- Number of chromosomes in the sample
- Length of the sequence
- Pi/sequence length 

Usage: "wisnoski_hw1.py sequence.fa"
"""

import sys
import math

# checks to ensure usage is correct
if len(sys.argv) < 2:
	print usage
else:	
	# reads filename from CLI to open fasta file for analysis
	in_file_name = sys.argv[1]
	print "Current file: " + in_file_name
	
	# open fasta file with filename provided in command line
	in_file = open(in_file_name, "r")
	
	# create empty list and dictionary to hold seqs
	in_seqs = {}
	current_seq = ""
	
	# cycle through lines in file and add to dictionary[ID] -> "sequence"
	for each in in_file:
		each = each.strip()
		if each.startswith(">"): # sequence ID
			current_seq = each
			in_seqs[current_seq] = "" # create dictionary key
		else:
			in_seqs[current_seq] += each # add sequence to dictionary 
	
	# calculate num chromosomes
	chrom_count = len(in_seqs.keys())
	
	seqs_by_position = {}
	position = 0 # index through each position to setup the list
	
	# initialize dictionary of nucleotides at each locus
	for each in in_seqs.values():
		position = 0 # index from the beginning again
		for base in each:
			seqs_by_position[position] = []
			position += 1
	
	# add values to the dictionary
	for each in in_seqs.values():
		position = 0 # index from the beginning again
		for base in each:
			seqs_by_position[position] += base
			position += 1
	
	# initialize key variables
	snp_count = 0
	pi = 0.0
	heterozygosities = []
	seq_lengths = []		
	seq_len = 0
	
	# Calculate Pi
	for locus in seqs_by_position.values(): # iterate over list of nucs at each position
		set_of_bases = set(locus) # to asses nucleotide diversity at each locus
		#print set_of_bases
		if len(set_of_bases) != 1: # i.e., there IS a SNP in this position
			#print locus
			nucleotides = set(['A','C','T','G'])
			if set_of_bases.issubset(nucleotides): # only contains A,C,T,G
				#print locus
				snp_count += 1 # update SNP counter
				pop_diversity = []
				seq_len += 1
				for base in set_of_bases: # for all bases present at this locus
					count = float(locus.count(base))
					pop_diversity.append(count/len(locus)) # add p1, p2, etc. to list
				squares = [p * p for p in pop_diversity] # square p1, p2, etc.				
				
				# calculate heterozygosity at this particular locus @ add to list
				h = (chrom_count/(chrom_count-1.0)) * (1 - sum(squares))
				heterozygosities.append(h)
				
			elif len(set_of_bases & nucleotides) != 1: # there is still a SNP here
				snp_count += 1
				pop_diversity = []
				indels = 0
				seq_len +=1
				for base in set_of_bases ^ nucleotides: # for indels present at this locus
					count = float(locus.count(base))
					indels += count # tracks how many indels at this locus
					#pop_diversity.append(count/len(locus)) # add p1, p2, etc. to list
				for base in set_of_bases & nucleotides: # for nucles present at this locus
					count = float(locus.count(base))
					pop_diversity.append(count/(len(locus)-indels)) # add p1, p2, etc. to list
				squares = [p * p for p in pop_diversity] # square p1, p2, etc.	
				
				# adjust heterozygosity for number of indels
				h = ((chrom_count-indels)/(chrom_count-1.0-indels)) * (1 - sum(squares))
				heterozygosities.append(h)
			elif len(set_of_bases & nucleotides) == 1:
				for each in set_of_bases & nucleotides:
					if each in "ACTG":
						if locus.count(each) > 1:
							seq_len += 1
		else:
			seq_len += 1	
					
	pi = sum(heterozygosities)

	# Calculate Wattersons Theta
	a1_i = []
	a2_i = []
	for i in range(chrom_count-1): # sums 1 / 1 - 11
		 a1_i.append(1.0/(i+1))
		 a2_i.append(1.0/((i+1)*(i+1)))
	a1 = sum(a1_i)
	a2 = sum(a2_i)
	theta_W = snp_count / a1

	# Calculate Tajima's D
	n = chrom_count
	d = pi - theta_W
	S = snp_count

	b1 = (n+1)/(3*(n-1.0))
	c1 = b1 - 1.0 / a1
	e1 = c1 / a1
	b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0))
	c2 = b2 - (n + 2.0)/(a1 * n) + a2/(a1*a1)
	e2 = c2 / (a1*a1 + a2)

	D = d / math.sqrt(e1*S + e2*S*(S-1))

	print "Pi: " + str(pi)
	print "Theta_w: " + str(theta_W)
	print "Tajima's D: " + str(D)

	
	# close file
	in_file.close()