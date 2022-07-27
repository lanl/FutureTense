#include "quasimodo.h"
#include "word.h"
#include <getopt.h>
#include <iostream>

using namespace std;

#define		DEFAULT_KMER_LEN			5
#define		DEFAULT_SEQUENCE_FILE		

Options::Options(int argc, char* argv[])
{
	bool print_usage = false;
	quit = false;

	kmer_len = DEFAULT_KMER_LEN;

	// Command line options:
	// [-i <input BLAST mapping file>] (default is stdin)
	// [-o <output file>] (default is stdout)
	// [-k <kmer length>] (default is 5)
	// [--ref <fasta file>]

	const char* options = "i:o:k:?h";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"ref", true, &config_opt, 1},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;

	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:

				if(config_opt == 1){ // ref

					reference_filename = optarg;
					break;
				}

				cerr << "Unknown command line flag!" << endl;
				
				quit = true;
				return;
				
			case 'i':
				input_mapping_filename = optarg;
				break;
			case 'o':
				output_filename = optarg;
				break;
			case 'k':
				kmer_len = abs( atoi(optarg) );
				break;
			case 'h':
			case '?':

				print_usage = true;
				break;
			default:
				cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
				
				quit = true;
				return;
		};
	}

	if(print_usage){

		cerr << "Usage for Quasimodo: " << QUASIMODO_VERSION << endl;
		cerr << "\t[-i <input SAM read-mapping file>] (default is stdin)" << endl;
		cerr << "\t[-o <output file>] (default is stdout)" << endl;
		cerr << "\t[-k <kmer length>] (default is " << DEFAULT_KMER_LEN << ")" << endl;
		cerr << "\t[--ref <fasta file>] (reference genome fasta file; single contig only)" << endl;

		quit = true;
		return;
	}

	if(optind < argc){

		cerr << "Unexpected additional command-line argument" << endl;
		quit = true;
		return;
	}

	if( (kmer_len < 1) || (kmer_len > MAX_WORD_LEN) ){

		cerr << "Please specify 1 <= kmer length <= " << MAX_WORD_LEN << endl;
		quit = true;
		return;
	}
}
