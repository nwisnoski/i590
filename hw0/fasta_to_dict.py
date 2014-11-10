#! /usr/bin/env python

usage = """
fasta_to_dict - Nathan Wisnoski
Reads a fasta file into a dictionary

Usage: fasta_to_dict.py sequence.fa
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
	
	print "Dictionary:"
	print in_seqs
	
	# close file
	in_file.close()