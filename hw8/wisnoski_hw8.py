#! /usr/bin/env python

usage = """
Homework 8 - Nathan Wisnoski
Calculates pi and 

usage: $ python wisnoski_hw8.py file.ms
where first line in file.ms startswith "./msdir" or "ms"
"""

import sys
import itertools

def parse_ms(in_file):
	parsed_ms = {}
	current_dataset = 0
	line_num = -3
	for line in in_file:
		line = line.strip()
		line_num += 1
		if line.startswith('ms'): # identify dataset characteristics
			line = line.split(' ')
			num_chromosomes = int(line[1])
			num_datasets = int(line[2])
		elif line.startswith('./'): # identify dataset characteristics
			line = line.split(' ')
			num_chromosomes = int(line[1])
			num_datasets = int(line[2])
		elif line == '//':
			current_dataset += 1 # new dataset
			parsed_ms[current_dataset] = []
		elif line_num > 0:
			dataset = parsed_ms[current_dataset]
			dataset.append(line)
			parsed_ms[current_dataset] = dataset
	
	return parsed_ms, num_datasets, num_chromosomes

def SNPs_by_position(snp_list):
	seq_length = len(snp_list[0])
	snp_dict = {}
	for position in range(seq_length):
		snp_dict[position] = []
	for chrom in snp_list:
		pos = 0
		for base in chrom:
			snp_dict[pos].append(base)
			pos += 1
	return snp_dict

def calc_pi(snp_dict, num_sites):
	
	h_values = []
	for position in snp_dict.keys():
		current_snp = snp_dict[position]
		num_chromosomes = float(len(current_snp))
		p_0 = float(current_snp.count('0'))/num_chromosomes
		p_1 = float(current_snp.count('1'))/num_chromosomes
		p0_squared = p_0 * p_0
		p1_squared = p_1 * p_1
		h = (num_chromosomes/(num_chromosomes-1.0))*(1.0-(p0_squared + p1_squared))
		h_values.append(h)
	
	pi = sum(h_values)

	return pi

def manage_calculations(parsed_data,num_datasets,num_chromosomes):
	current_dataset = 1
	pi_values = []
	while num_datasets - current_dataset > -1:
		# parse through each dataset to pull out relevant info
		this_set = parsed_data[current_dataset]
		num_sites = int(this_set[0][10:])
		positions = this_set[1][11:].split(' ')
		SNPs = this_set[2:]

		# calculate pi for each dataset
		SNP_dict = SNPs_by_position(SNPs)
		pi = calc_pi(SNP_dict, num_sites)
		
		pi_values.append(pi)
		
		current_dataset += 1
	
	# count number of times pi is less than or equal to 1
	count = 0
	for pi in pi_values:
		if pi < float(1):
		 	print pi
		 	count += 1
		elif float(str(pi)) == 1.0:
		 	print pi
		 	count += 1
	return count

if len(sys.argv) < 2:
	print usage
else:
	filename = sys.argv[1]
	in_file = open(filename, 'r')

	parsed_data, num_datasets, num_chromosomes = parse_ms(in_file)
	count = manage_calculations(parsed_data, num_datasets, num_chromosomes)

	print "The number of sites at which Pi is 1 or less is: " + str(count)
	print "so the resulting p-value is: " + str(float(count)/num_datasets)
	

	in_file.close()

