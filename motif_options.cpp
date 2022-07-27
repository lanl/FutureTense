#include "motif_options.h"
#include "motif.h"
#include <getopt.h>
#include <iostream>

using namespace std;

#define		DEFAULT_MOTIF_LEN		3
#define		DEFAULT_NUM_MOTIF		1
#define		DEFAULT_NUM_BATCH		1
#define		DEFAULT_NUM_EPOCH		1
#define		DEFAULT_NUM_FIT_TRIAL	1

#define		DEFAULT_REFERENCE_FILE	/*$HOME + */ "/ViralEvo/ref/MN908947.fna"
#define		DEFAULT_GENOMES_FILE	/*$HOME + */ "/ViralEvo/src/futuretense/sars_cov2_ncbi_june_29_2021.fna"
#define		DEFAULT_GENEALOGY_FILE	/*$HOME + */ "/ViralEvo/src/futuretense/parent_child.diff"
#define		DEFAULT_CODON_FILE		/*$HOME + */ "/ViralEvo/src/futuretense/sars_cov2_ncbi_june_29_2021.codon.gz"

void MotifOptions::load(int argc, char* argv[])
{
	bool print_usage = false;
	bool ignore_codon = false;

	quit = false;

	const char *home = getenv("HOME");

	if(home == NULL){
		throw __FILE__ ":MotifOptions::load: Unable to read the HOME environment variable";
	}
	
	reference_filename = string(home) + DEFAULT_REFERENCE_FILE;
	genomes_filename = string(home) + DEFAULT_GENOMES_FILE;
	genealogy_filename = string(home) + DEFAULT_GENEALOGY_FILE;
	codon_filename = string(home) + DEFAULT_CODON_FILE;

	motif_len = DEFAULT_MOTIF_LEN;
	num_motif = DEFAULT_NUM_MOTIF;
	num_batch = DEFAULT_NUM_BATCH;
	num_epoch = DEFAULT_NUM_EPOCH;
	num_trial = DEFAULT_NUM_FIT_TRIAL;
	cross_validation = TEST_ON_BRANCH_TIPS;
	distance_from_root = 0;
	seed = 0; // 0 --> use a time-based seed

	use_regions = false;

	// Command line options:
	// [--ref <reference genome filename>] (default is DEFAULT_REFERENCE_FILE)
	// [--genomes <genome sequence filename>] (default is DEFAULT_GENOMES_FILE)
	// [--genealogy <genealogy filename>] (default is DEFAULT_GENEALOGY_FILE)
	// [--codon <codon location filename>] (default is DEFAULT_CODON_FILE)
	// [--no-codon] (disable the use of codon location information)
	// [--regions] (enable region-specific motifs)
	// [--no-regions] (disable region-specific motifs)
	// [-l <motif length] (default is DEFAULT_MOTIF_LEN)
	// [-n <num motif>] (default is DEFAULT_NUM_MOTIF)
	// [-b <num batch>] (default is DEFAULT_NUM_BATCH)
	// [-e <num epoch>] (default is DEFAULT_NUM_EPOCH)
	// [-t <num fitting trials>] (default is DEFAULT_NUM_FIT_TRIAL)
	// [--eval <parameter file for objective function evaluation>]
	// [--cv.tip] (cross validate by testing on branch tips)
	// [--cv.root <distance from root>] (cross validate by testing on nodes >= distance from root)
	// [--seed <random number seed>] (0 uses time-based seed; default is 0)
	const char* options = "l:n:b:e:t:?h";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"cv.tip", false, &config_opt, 1},
		{"cv.root", true, &config_opt, 2},
		{"ref", true, &config_opt, 3},
		{"genomes", true, &config_opt, 4},
		{"genealogy", true, &config_opt, 5},
		{"codon", true, &config_opt, 6},
		{"no-codon", false, &config_opt, 7},
		{"regions", false, &config_opt, 8},
		{"no-regions", false, &config_opt, 9},
		{"seed", true, &config_opt, 10},
		{"eval", true, &config_opt, 11},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;

	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:

				if(config_opt == 1){ // cv.tip

					cross_validation = TEST_ON_BRANCH_TIPS;
					break;
				}

				if(config_opt == 2){ // cv.root

					cross_validation = TEST_BY_DISTANCE_FROM_ROOT;
					distance_from_root = abs( atoi(optarg) );
					break;
				}

				if(config_opt == 3){ // ref

					reference_filename = optarg;
					break;
				}

				if(config_opt == 4){ // genomes

					genomes_filename = optarg;
					break;
				}

				if(config_opt == 5){ // genealogy

					genealogy_filename = optarg;
					break;
				}

				if(config_opt == 6){ // codon

					codon_filename = optarg;
					break;
				}

				if(config_opt == 7){ // no-codon

					ignore_codon = true;
					break;
				}

				if(config_opt == 8){ // regions

					use_regions = true;
					break;
				}

				if(config_opt == 9){ // no-regions

					use_regions = false;
					break;
				}

				if(config_opt == 10){ // seed

					seed = abs( atoi(optarg) );
					break;
				}

				if(config_opt == 11){ // eval

					eval_filename.push_back(optarg);
					break;
				}

				cerr << "Unknown command line flag!" << endl;
				
				quit = true;
				return;
				
			case 'l':
				motif_len = abs( atoi(optarg) );
				break;
			case 'n':
				throw __FILE__ ":MotifOptions: Setting the number of motifs is not currently supported (-n)";
				num_motif = abs( atoi(optarg) );
				break;
			case 'b':
				num_batch = abs( atoi(optarg) );
				break;
			case 'e':
				num_epoch = abs( atoi(optarg) );
				break;
			case 't':
				num_trial = abs( atoi(optarg) );
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

		cerr << "Usage for Motif: " << MOTIF_VERSION << endl;

		cerr << "\t[--ref <reference genome filename>] (default is " << home << DEFAULT_REFERENCE_FILE << ")" << endl;
		cerr << "\t[--genomes <genome sequence filename>] (default is " << home << DEFAULT_GENOMES_FILE << ")" << endl;
		cerr << "\t[--genealogy <genealogy filename>] (default is " << home << DEFAULT_GENEALOGY_FILE << ")" << endl;
		cerr << "\t[--codon <codon location filename>] (default is " << home << DEFAULT_CODON_FILE << ")" << endl;
		cerr << "\t[--no-codon] (disable the use of codon location information)" << endl;
		cerr << "\t[--regions] (enable region-specific motifs)" << endl;
		cerr << "\t[--no-regions] (disable region-specific motifs; default)" << endl;

		cerr << "\t[-l <motif length] (default is " << DEFAULT_MOTIF_LEN << ")" << endl;
		//cerr << "\t[-n <num motif>] (default is " << DEFAULT_NUM_MOTIF << ")" << endl;
		cerr << "\t[-b <num batch>] (default is " << DEFAULT_NUM_BATCH << ")" << endl;
		cerr << "\t[-e <num epoch>] (default is " << DEFAULT_NUM_EPOCH << ")" << endl;
		cerr << "\t[-t <num fitting trials>] (default is " << DEFAULT_NUM_FIT_TRIAL << ")" << endl;
		cerr << "\t[--eval <parameter file for objective function evaluation> (may be repeated)" << endl;
		cerr << "\t[--cv.tip] (cross validate by testing on branch tips; default)" << endl;
		cerr << "\t[--cv.root <distance from root>] (cross validate by testing on nodes >= distance from root)" << endl;
		cerr << "\t[--seed <random number seed>] (0 uses time-based seed; default is 0)" << endl;
		
		quit = true;
		return;
	}

	if( (motif_len == 0) || (motif_len > MAX_MOTIF_LEN) ){

		cerr << "Please specify 1 <= motif length <= " << MAX_MOTIF_LEN << " (-l)"<< endl;

		quit = true;
		return;
	}

	if(motif_len%2 != 1){

		cerr << "Please specify an odd motif length (i.e., 1, 3, 5, ...; -l)"<< endl;

		quit = true;
		return;
	}

	if(num_motif == 0){

		cerr << "Please specify a motif redundancy of at least 1 (-n)" << endl;

		quit = true;
		return;
	}

	if(num_batch == 0){

		cerr << "Please specify at least one batch per epoch (-b)" << endl;

		quit = true;
		return;
	}

	if(num_epoch == 0){

		cerr << "Please specify at least one epoch (-e)" << endl;

		quit = true;
		return;
	}

	if(num_trial == 0){

		cerr << "Please specify at least one fitting trial (-t)" << endl;

		quit = true;
		return;
	}

	if( reference_filename.empty() ){

		cerr << "Please specify a reference genome file (--ref)" << endl;

		quit = true;
		return;
	}

	if( genomes_filename.empty() ){

		cerr << "Please specify a file of genome sequences (--genomes)" << endl;

		quit = true;
		return;
	}

	if( genealogy_filename.empty() ){

		cerr << "Please specify a genealogy file (--genealogy)" << endl;

		quit = true;
		return;
	}

	if( codon_filename.empty() ){

		cerr << "Please specify a file of codon locations (--codon)" << endl;

		quit = true;
		return;
	}

	if(ignore_codon){
		codon_filename = "";
	}
}