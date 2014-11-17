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

def calc_pi_theta(snp_dict, num_sites):
	h_values = {}

	for position in snp_dict.keys():
		current_snp = snp_dict[position]
		num_chromosomes = float(len(current_snp))
		p_0 = float(current_snp.count('0'))/num_chromosomes
		p_1 = float(current_snp.count('1'))/num_chromosomes
		p0_squared = p_0 * p_0
		p1_squared = p_1 * p_1
		h = (num_chromosomes/(num_chromosomes-1.0))*(1.0-(p0_squared + p1_squared))
		h_values[position] = h
		
	a_i = [1.0/(n+1.0) for n in range(num_sites)]
	a = sum(a_i)
	theta = num_sites / a
	pi = sum(h_values.values()) #/ num_sites

	return pi, theta

def calc_r_sqared(snp_dict, num_sites):
	positions = range(num_sites)
	snp_pairs = itertools.combinations(positions, 2)
	
	for pair in snp_pairs:
		snp1_pos, snp2_pos = pair[0], pair[1]
		snp1, snp2 = snp_dict[snp1_pos], snp_dict[snp2_pos]
		
		haplotypes = {}
		length = 0
		for each in snp1:
			length += 1

		for position in range(length):
			haplotypes[position] = [snp1[position], snp2[position]]
		unique_haps = {}
		denominator = 0

		allele_freqs1 = {}
		allele_freqs2 = {}

		for hap in haplotypes.values():
			hap_str = ''.join(hap)
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
	print '\n'

def manage_calculations(parsed_data,num_datasets,num_chromosomes):
	current_dataset = 1
	pi_values = []
	while num_datasets - current_dataset > -1:
		# parse through each dataset to pull out relevant info
		this_set = parsed_data[current_dataset]
		num_sites = int(this_set[0][10:])
		positions = this_set[1][11:].split(' ')
		SNPs = this_set[2:]

		SNP_dict = SNPs_by_position(SNPs)
		pi, theta = calc_pi_theta(SNP_dict, num_sites)
		
		pi_values.append(pi)

		# print "Dataset: " + str(current_dataset)
		# print "Pi = " + str(pi)
		# print "Watterson's Theta = " + str(theta)
		# calc_r_sqared(SNP_dict, num_sites)
		
		current_dataset += 1
	count = 0
	for each in pi_values:
		if float(each) <= float(1):
		 	print each
		 	count += 1
	return count



if len(sys.argv) < 2:
	print usage
else:
	filename = sys.argv[1]
	in_file = open(filename, 'r')

	parsed_data, num_datasets, num_chromosomes = parse_ms(in_file)
	count = manage_calculations(parsed_data, num_datasets, num_chromosomes)

	print "The number of sites at which Pi is < or == 1 is: " + str(count)
	print "so the resulting p-value is: " + str(count/num_datasets)

	in_file.close()

