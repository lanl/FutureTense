#!/usr/bin/env python3

# 1) Read a fasta file of SARS-CoV-2 genome sequences
# 2) Run blastx in the "blast two sequences" mode to align the nucleic acid genome sequences (as query)
#    to Wuhan protein sequences (as subject)
# 3) Parse the blastx output to produce the per-nucleotide genome position -> amino acid & codon position mapping
#    - Since the amino acid sequence can be recovered from the nucleic acid sequence, we only need to store the codon position
#      of each nucleotide:
#       - 0 = not coding sequence
#       - 1 = first codon position
#       - 2 = second codon position
#       - 3 = third codon position
#    - SARS-CoV-2 does not encode any amino acid sequence on the negative strand, so we don't need to store the sign of the reading frame
#    - Record the output as two column, tab-delimited: [genome accession][codon mapping]
#           accession1      00000000123123123123123123123123121231231230000012312312312312300000
#                                                             ^
#                                                             |
#                                                             Frame shift causes this nucleotide to be in both the
#                                                             third codon position of the upstream amino acid and in the
#                                                             first codon position of the downstream amino acid
import sys
import subprocess
import re

version = '0.1 Nov 22, 2021'

def main():
    
    subject_file = '/home/jgans/ViralEvo/ref/MN908947.faa' # <-- protein sequences from MN908947
    genome_file = './sars_cov2_ncbi_june_29_2021.fna'

    #####################################################################
    # Parse the genome fasta file
    print('Reading genome sequences', file=sys.stderr)

    fin = open(genome_file, 'r')

    if not fin:
        throw_error( "Unable to open {} for reading\n".format(genome_file) )
	
    accession = ''
    seq = ''
    count = 0

    for line in fin:

        # Skip the header
        if '>' in line:
			
            # Extract the sequence accession from the fasta header
            match = re.search('>\s*(\w+)\s', line)

            if match:

                if seq != '':
                    extract_codons(accession, seq, subject_file)
                    count += 1

                accession = match.group(1)

                # 'reference' is a special synonym for SARS-CoV-2 Wuhan reference genome
                if accession == 'MN908947':
                    accession = 'reference'
            
            else:
                throw_error( "Unable to extract accession from {}\n".format(line) )
            
            seq = ''

        else:
            seq += line.strip()

    fin.close()

    if seq != '':
        extract_codons(accession, seq, subject_file)
        count += 1

    print('Wrote codon mapping for {} sequences'.format(count), file=sys.stderr)

# Write the codon map to stdout
def extract_codons(m_accession, m_seq, m_subject_file):

    # Run blastx in blast-two-sequences mode
    command = ['/home/jgans/ncbi-blast-2.12.0+-src/bin/blastx', 
        '-query', '-', 
        '-subject', m_subject_file, 
        '-seg', 'no',
        '-strand', 'plus',
        '-evalue', str(0.1)]

    p = subprocess.Popen(command,
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE)

    # Pipe the input genome sequence vi stdin
    # communicate returns (stdout, stderr)
    blast_out = p.communicate( input = bytearray('>{}\n{}'.format(m_accession, m_seq), 'utf-8') )[0].decode("utf-8")

    codon = [0]*len(m_seq)

    # Split the blastx output by lines
    for line in blast_out.split('\n'):
        
        m = re.search('^Query\s+(\d+)\s+\w+\s+(\d+)$', line)

        if m:
            
            # Convert from 1's based BLAST output to zero based indicies
            start = int(m.group(1)) - 1
            stop = int(m.group(2)) - 1

            loc = 0

            for i in range(start, stop + 1):
                
                codon[i] = loc + 1 # Store codon position as 1's based numbers

                loc = (loc + 1)%3

    print('{}\t{}'.format(m_accession, ''.join( str(x) for x in codon)) )

def throw_error(m_error):

    print(m_error, file=sys.stderr)
    sys.exit(1)


main()