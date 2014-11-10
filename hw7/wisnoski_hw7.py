#! /usr/bin/env python

usage = """
Homework 7 - Nathan Wisnoski
Calculates pi, theta, and r^2 from each dataset.

usage: $ python wisnoski_hw7.py file.ms
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
	positions = range(num_sites)
	snp_pairs = itertools.combinations(positions, 2)
	print snp_pairs
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
	pi = sum(h_values.values()) / num_sites


	return pi, theta


def calculations(parsed_data,num_datasets,num_chromosomes):
	current_dataset = 1
	while num_datasets - current_dataset > -1:
		# parse through each dataset to pull out relevant info
		this_set = parsed_data[current_dataset]
		num_sites = int(this_set[0][10:])
		positions = this_set[1][11:].split(' ')
		SNPs = this_set[2:]

		SNP_dict = SNPs_by_position(SNPs)
		pi, theta = calc_pi_theta(SNP_dict, num_sites)
		
		print "Dataset: " + str(current_dataset)
		print "Pi = " + str(pi)
		print "Watterson's Theta = " + str(theta) + '\n'

		current_dataset += 1


if len(sys.argv) < 2:
	print usage
else:
	filename = sys.argv[1]
	in_file = open(filename, 'r')

	parsed_data, num_datasets, num_chromosomes = parse_ms(in_file)
	calculations(parsed_data, num_datasets, num_chromosomes)

	in_file.close()

