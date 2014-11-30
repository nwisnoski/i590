#! /usr/bin/env python

usage = """
Homework 8 - Nathan Wisnoski
Calculates pi for .ms file and returns
pi for each dataset and the p-value for pi=1

usage: $ python wisnoski_hw8.py file.ms
where first line in file.ms startswith "./msdir" or "ms"
"""

import sys
import itertools
import math

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

def calc_pi_theta_D(snp_dict, num_sites):
	# calculate Pi
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

	n = num_chromosomes
	# Calculate Wattersons Theta
	a1_i = []
	a2_i = []
	for i in range(int(n)-1): # sums 1 / 1 - 11
		 a1_i.append(1.0/(i+1))
		 a2_i.append(1.0/((i+1)*(i+1)))
	a1 = sum(a1_i)
	a2 = sum(a2_i)
	theta_W = num_sites / a1

	# calculate tajima's d
	d = pi - theta_W
	S = num_sites

	b1 = (n+1)/(3*(n-1.0))
	c1 = b1 - 1.0 / a1
	e1 = c1 / a1
	b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0))
	c2 = b2 - (n + 2.0)/(a1 * n) + a2/(a1*a1)
	e2 = c2 / (a1*a1 + a2)

	D = d / math.sqrt(e1*S + e2*S*(S-1))

	return D

def manage_calculations(parsed_data,num_datasets,num_chromosomes):
	current_dataset = 1
	pi_values = []
	D_values = []
	while num_datasets - current_dataset > -1:
		# parse through each dataset to pull out relevant info
		this_set = parsed_data[current_dataset]
		num_sites = int(this_set[0][10:])
		positions = this_set[1][11:].split(' ')
		SNPs = this_set[2:]

		# calculate pi for each dataset
		SNP_dict = SNPs_by_position(SNPs)
		D = calc_pi_theta_D(SNP_dict, num_sites)
		
		# pi_values.append(pi)
		D_values.append(D)
		
		current_dataset += 1
	
	# count number of times pi is less than or equal to 1
	count = 0
	for D in D_values:
		print D
		if D < -1.64462926178:
		 	count += 1
		elif float(str(D)) == -1.64462926178:
		 	count += 1
	return count

if len(sys.argv) < 2:
	print usage
else:
	filename = sys.argv[1]
	in_file = open(filename, 'r')

	parsed_data, num_datasets, num_chromosomes = parse_ms(in_file)
	count = manage_calculations(parsed_data, num_datasets, num_chromosomes)

	print "The number of sites at which D is -1.64462926178 or less is: " + str(count)
	print "so the resulting p-value is: " + str(float(count)/num_datasets)
	

	in_file.close()

