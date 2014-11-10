#! /usr/bin/env python

usage = """
count_bases.py - Nathan Wisnoski
Reads a fasta file into a dictionary and displays
a count of each nucleotide base

Usage: count_bases.py sequence.fa
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
	
	# cycle through lines in file
	for each in in_file:
		each = each.strip()
		if each.startswith(">"): # sequence ID
			current_seq = each
			in_seqs[current_seq] = "" # create dictionary key
		else:
			in_seqs[current_seq] += each # add sequence to dictionary 
	
	# Toggle print sequences on/off to check dictionary
	#print "Sequence dictionary:"
	#print in_seqs
	
	# create dictionary to count base composition & initialize to 0
	nucleotides = {}
	for base in "ACTG":
		nucleotides[base] = 0
	
	# count bases
	for seq in in_seqs.values(): # loops through sequences only
		for base in list(set(seq)):
			count = seq.count(base)
			nucleotides[base] += count
	
	for base in "ACTG":
		print base + ": " + str(nucleotides[base])
	
	# close file
	in_file.close()