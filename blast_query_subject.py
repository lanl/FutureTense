#!/usr/bin/env python3

# 1) Read a *.sif file of parent child relationships
# 2) Read a fasta file of genome sequences (whose accessions are listed in the *.sif file)
# 3) Read the "reference" genome
# 4) Run blastn in the "blast two sequences" mode to align the child (as query) to the parent (as subject)

import sys
import subprocess
import re

version = '0.1 Aug 25, 2021'

def main():

	temp_subject_file = '_subject.fna'
	temp_query_file = '_query.fna'

	input_sif_file = './sars_cov2_ncbi_june_29_2021.sif'
	input_reference_fasta = '/Users/jgans/ViralEvo/ref/MN908947.fna'
	input_fasta = './sars_cov2_ncbi_june_29_2021.fna'

	seq = dict() # accession -> sequence

	#####################################################################
	# Parse the reference fasta file
	print('Reading reference genome', file=sys.stderr)

	fin = open(input_reference_fasta, 'r')

	if not fin:

		sys.stderr.write( "Unable to open {} for reading\n".format(input_reference_fasta) )
		sys.exit(1)
	
	# Please note that 'reference' is a special accession name
	seq['reference'] = ''

	for line in fin:

		# Skip the header
		if '>' in line:
			continue
		else:

			# Please note that 'reference' is a special accession name
			seq['reference'] += line.strip()

	fin.close()

	#####################################################################
	# Parse the genome fasta file
	print('Reading genome sequences', file=sys.stderr)

	fin = open(input_fasta, 'r')

	if not fin:

		sys.stderr.write( "Unable to open {} for reading\n".format(input_fasta) )
		sys.exit(1)
	
	accession = ''

	for line in fin:

		# Skip the header
		if '>' in line:
			
			# Extract the sequence accession from the fasta header
			match = re.search('>\s*(\w+)\s', line)

			if match:

				accession = match.group(1)
			else:
				sys.stderr.write( "Unable to extract accession from {}\n".format(line) )
				sys.exit(1)

			seq[accession] = ''
		else:
			seq[accession] += line.strip()

	fin.close()

	#####################################################################
	# Parse the sif file
	print('Reading parent-child relationships', file=sys.stderr)
	
	fin = open(input_sif_file, 'r')

	if not fin:

		sys.stderr.write( "Unable to open {} for reading\n".format(input_sif_file) )
		sys.exit(1)
	
	for line in fin:

		# Skip any comments
		if '#' in line:
			line = line[:line.index('#')]

		data = line.strip().split('\t')

		if len(data) != 3:
			continue

		if data[1] != 'parent_of':

			sys.stderr.write( "Invalid relationship in {}\n".format(input_sif_file) )
			sys.exit(1)

		# The parent (subject) accession is data[0]

		if data[0] not in seq:
			
			sys.stderr.write( "Unable to lookup parent (subject) accession: \'{}\'\n".format(data[0]) )
			sys.exit(1)

		# The child (query) accession is data[2]
		if data[2] not in seq:
			
			sys.stderr.write( "Unable to lookup child (query) accession: \'{}\'\n".format(data[2]) )
			sys.exit(1)

		# Write the parent (subject)
		fout = open(temp_subject_file, 'w')

		if not fout:

			sys.stderr.write( "Unable to open {} for writing\n".format(temp_subject_file) )
			sys.exit(1)
		
		fout.write('>{}\n'.format(data[0]))
		fout.write('{}\n'.format(seq[ data[0] ]))

		fout.close()

		# Write the child (query)
		fout = open(temp_query_file, 'w')

		if not fout:

			sys.stderr.write( "Unable to open {} for writing\n".format(temp_query_file) )
			sys.exit(1)
		
		fout.write('>{}\n'.format(data[2]))
		fout.write('{}\n'.format(seq[ data[2] ]))

		fout.close()

		# Run blastn in blast-two-sequences mode
		command = ['/Users/jgans/src/tntblast-2.2/ncbi-blast-2.11.0+-src/bin/blastn', 
			'-query', temp_query_file, 
			'-subject', temp_subject_file, 
			'-dust', 'no']

		ret = subprocess.run(command, capture_output=True)

		sys.stdout.write( ret.stdout.decode("utf-8") )

	fin.close()
	
main()