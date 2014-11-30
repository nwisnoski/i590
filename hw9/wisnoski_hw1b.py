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
	

	
	
		
		#print in_seqs[seq_id]
	
	
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
					
	# here, heterozygosities should be a list of h at each location in our sequence
	#print heterozygosities
	pi = sum(heterozygosities)
	
	print "Number of SNPs: " + str(snp_count)
	print "Number of chromosomes: " + str(chrom_count)
	print "Sequence length: " + str(seq_len)
	print "Pi/sequence: " + str(pi/seq_len)
	
	# close file
	in_file.close()