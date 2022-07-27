#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <list>

#include <stdlib.h>
#include <mpi.h>

#include "parse_fasta.h"
#include "parse_codon.h"
#include "genome_difference.h"
#include "keys.h"
#include "shuffle.h"
#include "motif.h"
#include "word.h"
#include "map.h"
#include "mpi_util.h"
#include "motif_options.h"
#include "tensor.h"
#include "index.h"
#include "split.h"

using namespace std;

// FYI regarding std::numeric_limits:
//
// 		std::numeric_limits<float>::min() = 1.17549e-38 <-- Smallest positive, non-zero number
// 		std::numeric_limits<float>::max() = 3.40282e+38 <-- largest positive number
// 		std::numeric_limits<float>::lowest() = -3.40282e+38 <-- Negative number with the largest absolute value

// Global variables for MPI
int mpi_numtasks;
int mpi_rank;

#define		IS_ROOT		(mpi_rank == 0)

#define		EPS			1.0e-7f

// All of the tensor dimensional parameters in one place!
struct TensorDim
{
	unsigned int num_codon_loc;
	unsigned int num_region;
	unsigned int num_mutation;
};

void load_genealogy(MULTIMAP<string /*parent*/, string /*child*/> &m_family, MAP<string /*child*/, GenomeDifference> &m_diff,
	const string &m_filename);
Tensor<float> find_motifs(const TensorDim &m_var, const deque< pair<string, string> > &m_train, const MAP<string, Sequence> &m_db, 
	const MAP<string, GenomeDifference> &m_diff, const vector< pair<unsigned int, unsigned char> > &m_regions,
	const MotifOptions &m_opt, const gsl_rng *m_rand_ptr);
float score(const vector< deque<WordIndex> > &m_mut, const Tensor<float> &m_buffer, const unsigned int &m_kmer_len);
deque< pair<float /*child accession*/, string /*child auroc*/> >  
	predict_mutations(const deque< Tensor<float> > &m_param, const TensorDim &m_var,
	const deque< pair<string, string> > &m_test, const MAP<string, Sequence> &m_db, 
	const MAP<string, GenomeDifference> &m_diff, const vector< pair<unsigned int, unsigned char> > &m_regions,
	const unsigned int &m_motif_len);
vector<bool> mask_secondary_structure(const Sequence &m_seq, const unsigned int &m_kmer_len,
	const unsigned int &m_window);
size_t distance_to_root(const string &m_node, const MAP<string, string> &m_child_to_parent);
unsigned char get_region(const vector< pair<unsigned int /*zero-based start*/, unsigned char /*region id*/> > &m_regions, 
	const unsigned int &m_loc);
const pair<Word, bool> extract_centered_kmer(const unsigned int &m_center_loc, const unsigned int &m_motif_len, 
	const Sequence &m_seq);
Tensor<float> load_param(const string &m_filename, const TensorDim &m_dim, const size_t &m_num_word);

template<class CONTAINER> double average_first(const CONTAINER &m_begin, const CONTAINER &m_end)
{
	size_t norm = 0;
	double ret = 0.0;

	for(CONTAINER i = m_begin;i != m_end;++i){

		++norm;
		ret += i->first;
	}

	if(norm == 0){
		return 0.0;
	}

	return ret/norm;
}

int main(int argc, char* argv[])
{
	try{

		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

		MotifOptions opt;

		if(IS_ROOT){
			opt.load(argc, argv);
		}

		broadcast(opt, mpi_rank, 0);

		 if(opt.quit){

			MPI_Finalize();
			return EXIT_SUCCESS;
		}

		double profile = MPI_Wtime();

		if(IS_ROOT){

			cerr << "Version " << MOTIF_VERSION << endl;
			cerr << "Running with " << mpi_numtasks << " MPI ranks" << endl;
			cerr << "Reference genome = " << opt.reference_filename << endl;
			cerr << "Genome sequences = " << opt.genomes_filename << endl;
			cerr << "Genealogy = " << opt.genealogy_filename << endl;

			if( !opt.codon_filename.empty() ){
				cerr << "Per nucleotide codon positions = " << opt.codon_filename << endl;
			}

			cerr << "Motif length = " << opt.motif_len << endl;
			cerr << "Number of batches = " << opt.num_batch << endl;
			cerr << "Number of epochs = " << opt.num_epoch << endl;
			cerr << "Number of fitting trials = " << opt.num_trial << endl;

			switch(opt.cross_validation){
				case TEST_ON_BRANCH_TIPS:
					cerr << "Cross validation by testing on branch tip children" << endl;
					break;
				case TEST_BY_DISTANCE_FROM_ROOT:
					cerr << "Cross validation by testing on children that are >= " 
						<< opt.distance_from_root << " levels from the tree root" << endl;
					break;
				default:
					cerr << "Unknown cross validation scheme!" << endl;
					throw __FILE__ ":main: Unknown cross validation scheme";
			};
		}

		// Initialize a random number generator
		gsl_rng *rand_gen = gsl_rng_alloc(gsl_rng_ranlux389);

		if(opt.seed == 0){

			// A seed of zero triggers a time-based seed
			opt.seed = time(NULL);
		}
		
		// Make sure all ranks use a unique random number seed
		broadcast(opt.seed, mpi_rank, 0);

		opt.seed += mpi_rank;

		gsl_rng_set(rand_gen, opt.seed);

		if(IS_ROOT){
			cerr << "Random number seed = " << opt.seed << endl;
		}
		
		// All of the genome sequences
		MAP<string /*accession*/, Sequence> db;

		// To avoid overloading simple NFS storage devices, only the root reads from disk
		if(IS_ROOT){

			// All of the genome sequences
			MAP<string /*accession*/, Sequence> reference_db;

			// Read the single reference genome 
			load_fasta(reference_db, opt.reference_filename);

			if(reference_db.size() != 1){
				throw __FILE__ ":main: Failed to read the expected single reference genome";
			}

			cerr << "The genome " << reference_db.begin()->first << " is the reference and contains "
				<< reference_db.begin()->second.size() << " nucleotides" << endl;

			// Add the reference genome accession to database with the special accession `reference`
			db["reference"] = reference_db.begin()->second;

			// Count the nucleotide frequency in the reference genome
			Tensor<float> composition(NUM_BASE);

			for(Sequence::const_iterator i = reference_db.begin()->second.begin();i != reference_db.begin()->second.end();++i){

				switch(i->base()){
					case SequenceBase::A:
						++composition(BASE_A);
						break;
					case SequenceBase::T:
						++composition(BASE_T);
						break;
					case SequenceBase::G:
						++composition(BASE_G);
						break;
					case SequenceBase::C:
						++composition(BASE_C);
						break;
				};
			}

			cerr << "Reference genome base composition: |A| = " << composition(BASE_A) 
				<< ", |C| = " << composition(BASE_C) 
				<< ", |G| = " << composition(BASE_G) 
				<< ", |T| = " << composition(BASE_T) << endl;

			cerr << "Reading sequences ... ";
		
			reference_db.clear(); // No longer used

			// Read the rest of the genome sequences
			load_fasta(db, opt.genomes_filename);

			cerr << "done." << endl;
		
			// The reference genome is MN908947, remove this accession if it also appears in the genome database
			if(db.erase("MN908947") > 0){

				if(IS_ROOT){
					cerr << "Removed MN908947 from the database (redundant with \"reference\")" << endl;
				}
			}

			cerr << "The genome database contains a total of " << db.size() << " sequences" << endl;

			if( !opt.codon_filename.empty() ){

				cerr << "Reading per-nucleotide codon mapping ... ";
				
				insert_codon_loc(opt.codon_filename, db);

				cerr << "done." << endl;
			}
			else{
				cerr << "No per-nucleotide codon mapping provided" << endl;
			}
		}

		// Share the sequence database to all nodes
		broadcast(db, mpi_rank, 0);
		
		TensorDim var;

		// There are four codon positions: 0 == non-coding, 1 == first codon position, 2 == second codon position and
		// 3 == third codon position
		var.num_codon_loc = opt.codon_filename.empty() ? 1 : 4;
		var.num_mutation = NUM_MUTATION;

		// The difference between the child and its parent. The parent is the reference frame
		MAP<string /*child accession*/, GenomeDifference> diff;
		MULTIMAP<string /*parent accession*/, string /*child accession*/> family;

		if(IS_ROOT){

			load_genealogy(family, diff, opt.genealogy_filename);

			cerr << "Found " << family.size() << " parent-child relationships (skipping redundant strains)" << endl;
		}

		// Share the genealogy data to all nodes
		broadcast(family, mpi_rank, 0);
		broadcast(diff, mpi_rank, 0);

		// The keys of the family map are all of the nodes that have at least one child
		const deque<string> parent_nodes = keys(family);

		if(IS_ROOT){
			cerr << "Found " << parent_nodes.size() << " nodes with one or more children" << endl;
		}

		// Invert the family multimap of parent->children to a *map* of child->parent (each child has one parent).
		// This is needed for computing the distance distribution of branch lengths.
		MAP<string /*child*/, string /*parent*/> inv_family;

		for(MULTIMAP<string, string>::const_iterator i = family.begin();i != family.end();++i){
			inv_family[i->second] = i->first;
		}

		// Compute the distribution of tip-to-root distances
		MAP<size_t /*distance to root*/, size_t /*number of nodes*/> distance_hist;

		for(MAP<string, string>::const_iterator i = inv_family.begin();i != inv_family.end();++i){
			++distance_hist[distance_to_root(i->first, inv_family)];
		}

		if(IS_ROOT){

			cerr << "Histogram: [distance to root] [number of nodes]" << endl;

			for(size_t d = 1;true;++d){

				MAP<size_t, size_t>::const_iterator iter = distance_hist.find(d);

				if( iter == distance_hist.end() ){
					break;
				}

				cerr << d << '\t' << iter->second << endl;
			}
		}

		// Define one or more genome regions. Only defining zero regions will cause an error
		// Region definition needs to be via input file, rather than hard coded!
		MAP<string, unsigned char> region_id;
		vector< pair<unsigned int /*zero-based start*/, unsigned char /*region id*/> > regions;
		
		var.num_region = 1;

		if(opt.use_regions){

			region_id["5' UTR"] = 0;
			region_id["orf1ab"] = 1;
			region_id["S"] = 2;
			region_id["ORF3a"] = 3;
			region_id["E"] = 4;
			region_id["M"] = 5;
			region_id["ORF6"] = 6;
			region_id["ORF7a"] = 7;
			region_id["ORF8"] = 8;
			region_id["N"] = 9;
			region_id["ORF10"] = 10;
			region_id["3' UTR"] = region_id["5' UTR"];

			if(IS_ROOT){

				cerr << "Region definitions:" << endl;

				for(MAP<string, unsigned char>::const_iterator i = region_id.begin();i != region_id.end();++i){
					cerr << '\t' << i->first << " -> " << int(i->second) << endl;
				}
			}

			// Region definitions from the MN908947 reference sequence GBK file.
			regions.push_back( make_pair(0, region_id["5' UTR"] ) );
			regions.push_back( make_pair(265, region_id["orf1ab"] ) );
			regions.push_back( make_pair(21562, region_id["S"] ) );
			regions.push_back( make_pair(25392, region_id["ORF3a"] ) );
			regions.push_back( make_pair(26244, region_id["E"] ) );
			regions.push_back( make_pair(26522, region_id["M"] ) );
			regions.push_back( make_pair(27201, region_id["ORF6"] ) );
			regions.push_back( make_pair(27393, region_id["ORF7a"] ) );
			regions.push_back( make_pair(27893, region_id["ORF8"] ) );
			regions.push_back( make_pair(28273, region_id["N"] ) );
			regions.push_back( make_pair(29557, region_id["ORF10"] ) );
			regions.push_back( make_pair(29674, region_id["3' UTR"] ) );

			// Must be sorted in ascending order to be compatible with lower_bound
			sort( regions.begin(), regions.end() );

			for(MAP<string, unsigned char>::const_iterator i = region_id.begin();i != region_id.end();++i){
				var.num_region = max(var.num_region, i->second + 1U);
			}
		}

		if(IS_ROOT){
			cerr << "Found a total of " << var.num_region << " regions" << endl;
		}

		// What is an appropriate null model for mutation prediction? A survey of GenBank SARS-CoV-2 from August 2021 data shows:
		//		Odds of substituting a parent base to any child base:
		//			p(A) = 0.116231786
		//			p(T) = 0.166989398
		//			p(G) = 0.230628843
		//			p(C) = 0.486149973
		// 		
		//			Conditional probability of a child base given a parent base: P(Child|Parent)
		//				A			T			G			C
		//			A	-			0.151228655	0.706747189	0.142024157
		//			T	0.123264241	-			0.150514567	0.726221191
		//			G	0.291597574	0.629142965	-			0.079259461
		//			C	0.044292201	0.944375846	0.011331953	-		
		//
		//			The rows are the parent base and the columns are the child base. Rows sum to 1.
		//
		//	From the above statistics, a potential null model for predicting evolution is to randomly select a "C" in the input
		//	genome and mutate it into a "T".
		//
		// Task #1: Find motifs in the nodes with one or more children and test the prediction accuracy of these motifs
		// on the leaf nodes (that do *not* have any children).

		// "Unpack" the training and testing parent-child pairs
		deque< pair<string /*parent*/, string /*child*/> > train;
		deque< pair<string /*parent*/, string /*child*/> > test;

		size_t num_skipped = 0;

		for(MULTIMAP<string /*parent*/, string /*child*/>::const_iterator i = family.begin();i != family.end();++i){

			// Exclude parent-child pairs in which the child differs from the parent by a non-ATGC base
			MAP<string, GenomeDifference>::const_iterator diff_iter = diff.find(i->second);

			if( diff_iter == diff.end() ){
				throw __FILE__ ":main: Unable to lookup child accession in diff";
			}

			if( !diff_iter->second.is_ATGC() ){

				++num_skipped;
				continue;
			}

			if(opt.cross_validation == TEST_ON_BRANCH_TIPS){

				// Is this a "terminal" parent-child pair? In this context, "terminal" means that the child accession is *not*
				// also a parent
				if( family.find(i->second) == family.end() ){

					// Terminal pair for testing
					test.push_back( make_pair(i->first, i->second) );
				}
				else{
					train.push_back( make_pair(i->first, i->second) );
				}
			}
			else if(opt.cross_validation == TEST_BY_DISTANCE_FROM_ROOT){

				if(distance_to_root(i->second, inv_family) >= opt.distance_from_root){
					test.push_back( make_pair(i->first, i->second) );
				}
				else{
					train.push_back( make_pair(i->first, i->second) );
				}
			}
			else{
				throw __FILE__ ":main: Unknown cross validation scheme";
			}
		}

		// Make sure that all ranks agree on the test and training set membership
		broadcast(test, mpi_rank, 0);
		broadcast(train, mpi_rank, 0);

		const size_t num_train = train.size();

		if(IS_ROOT){

			cerr << "Found " << num_train << " training parent-child pairs" << endl;
			cerr << "Found " << test.size() << " testing parent-child pairs" << endl;
			cerr << num_skipped << " pairs were skipped due to non-ATGC bases" << endl;
		}

		const Word num_word = 1 << (2*opt.motif_len) /* 4^L = 2^(2*L)*/;
		const Word central_base_offset = opt.motif_len/2;

		// What are the observed, parent to child substitution frequencies as a function of
		// parent codon position (i.e., 0, 1, 2, 3). Note that a codon position of 0 indicates that
		// the given nucleotide is not located in coding sequence. All SARS-CoV-2 proteins are encoded on
		// the plus strand RNA.
		Tensor<float> mutation_freq(var.num_codon_loc, var.num_region, 
			NUM_MUTATION /*child mutation*/, num_word /*parent sequence*/);
		
		//for(deque<string>::const_iterator p = parent_nodes.begin();p != parent_nodes.end();++p){
		for(deque< pair<string /*parent*/, string /*child*/> >::const_iterator iter = train.begin();iter != train.end();++iter){

			const string &parent = iter->first;
			const string &child = iter->second;

			MAP<string, Sequence>::const_iterator parent_iter = db.find(parent);

			if( parent_iter == db.end() ){
				throw __FILE__ ":main: Unable to look up parent sequence";
			}

			MAP<string, GenomeDifference>::const_iterator diff_iter = diff.find(child);

			if( diff_iter == diff.end() ){
				throw __FILE__ ":main: Unable to look up child GenomeDifference";
			}

			// Enumerate all of the substitutions that differentiate the parent from the child
			for(deque< pair<unsigned int, char> >::const_iterator i = diff_iter->second.get_substitution().begin();
				i != diff_iter->second.get_substitution().end();++i){

				// This gets a bit messy since sequences are stored with a 4 bit encoding and we need to map to
				// a 2-bit encoding for mutation_freq (and the genome differences are stored using ASCII 
				// characters!). Extract the kmer of length opt.motif_len centered at i->first.
				const pair<Word, bool> parent_word = extract_centered_kmer(i->first, opt.motif_len, parent_iter->second);

				if(parent_word.second == false){

					// Could not extract a valid word (due to degenerate base or running past the edge of the sequence)
					continue;
				}

				unsigned char child_mutation = NO_MUTATION;

				switch(i->second){
					case 'A':
						child_mutation = MUTATION_SUB_A;
						break;
					case 'T':
						child_mutation = MUTATION_SUB_T;
						break;
					case 'G':
						child_mutation = MUTATION_SUB_G;
						break;
					case 'C':
						child_mutation = MUTATION_SUB_C;
						break;
				};

				if(child_mutation == NO_MUTATION){
					continue;
				}

				++mutation_freq( parent_iter->second[i->first].codon(), get_region(regions, i->first), 
					child_mutation, parent_word.first);
			}

			#ifdef MUTATION_INSERTION
			for(deque< pair<unsigned int, string> >::const_iterator i = diff_iter->second.get_insertion().begin();
				i != diff_iter->second.get_insertion().end();++i){

				// This gets a bit messy since sequences are stored with a 4 bit encoding and we need to map to
				// a 2-bit encoding for mutation_freq (and the genome differences are stored using ASCII 
				// characters!).Extract the kmer of length opt.motif_len centered at i->first.
				const pair<Word, bool> parent_word = extract_centered_kmer(i->first, opt.motif_len, parent_iter->second);

				if(parent_word.second == false){

					// Could not extract a valid word (due to degenerate base or running past the edge of the sequence)
					continue;
				}

				++mutation_freq( parent_iter->second[i->first].codon(), get_region(regions, i->first), 
					MUTATION_INSERTION, parent_word.first);
			}
			#endif // MUTATION_INSERTION

			#ifdef MUTATION_DELETION
			for(deque<unsigned int>::const_iterator i = diff_iter->second.get_deletion().begin();
				i != diff_iter->second.get_deletion().end();++i){

				// This gets a bit messy since sequences are stored with a 4 bit encoding and we need to map to
				// a 2-bit encoding for mutation_freq (and the genome differences are stored using ASCII 
				// characters!). Extract the kmer of length opt.motif_len centered at i->first.
				const pair<Word, bool> parent_word = extract_centered_kmer(*i, opt.motif_len, parent_iter->second);

				if(parent_word.second == false){

					// Could not extract a valid word (due to degenerate base or running past the edge of the sequence)
					continue;
				}

				++mutation_freq( parent_iter->second[*i].codon(), get_region(regions, *i), 
					MUTATION_DELETION, parent_word.first);
			}
			#endif // MUTATION_DELETION
		}

		if(IS_ROOT){

			for(unsigned int c = 0;c < var.num_codon_loc;++c){

				for(unsigned int r = 0;r < var.num_region;++r){

					cerr << "Parent -> child base mutation frequency at codon position " << c << " and region " << r << endl;

					for(int i = 0;i < NUM_MUTATION;++i){
						cerr << '\t' << mutation_name[i];
					}

					cerr << endl;

					for(Word w = 0;w < num_word;++w){

						cerr << word_to_string(w, opt.motif_len);

						for(int j = 0;j < NUM_MUTATION;++j){
							cerr << '\t' << mutation_freq(c, r, j, w);
						}

						cerr << endl;
					}
				}
			}
		}

		// For a baseline (i.e., null) model of evolution prediction for SARS-CoV-2, we observe
		// that the most commonly observed mutation is a C->T substitution. Since the SARS-CoV-2
		// genome has a low, 38% G+C content, the 'C' is a rare base. A reasonable baseline model is
		// to precict all 'C' bases in the parent genome as being equally likely to transition to a 'T'
		// (and score all other parent bases and mutations with a lower score).
		deque< Tensor<float> > baseline( 1 /*ensemble size is 1*/,
			Tensor<float>(var.num_codon_loc, 
				var.num_region, 
				var.num_mutation,
				num_word) );
		
		// Populate all tensor elements corresponding to kmers with a central C nucleotide
		baseline[0] = 0.0;

		// See Simmonds, "Rampant C->U Hypermutation in the Genomes of SARS-CoV-2 and Other Coronaviruses: 
		// Causes and Consequences for Their Short- and Long-Term Evolutionary Trajectories", mSphere, 2020, vol. 5, no. 3
		for(unsigned int c = 0;c < var.num_codon_loc;++c){

			for(unsigned int r = 0;r < var.num_region;++r){

				for(Word w = 0;w < num_word;++w){

					// Is the central base in this word a 'C'?
					if( ( (w >> 2*central_base_offset) & 3 ) == BASE_C){

						// Since the predominantly observed mutation in SARS-CoV-2 is C->T,
						// assign all table words that contain a C in the central position a
						// mutation rank of 1 for mutating to a T.
						baseline[0](c, r, MUTATION_SUB_T, w) = 1.0;
					}
				}
			}
		}

		deque< pair<float, string> > baseline_scores = predict_mutations(baseline, var, test, db, diff, regions,
			opt.motif_len);

		float ave_baseline_score = average_first( baseline_scores.begin(), baseline_scores.end() );

		baseline_scores.clear();

		if(IS_ROOT){
			cerr << "Baseline model (all C->T mutations are equally likely, all other mutations equally unlikely) score = " 
				<< ave_baseline_score << endl;
		}

		// Used the observed distribution as baseline
		for(unsigned int c = 0;c < var.num_codon_loc;++c){
			for(unsigned int r = 0;r < var.num_region;++r){
				for(Word w = 0;w < num_word;++w){			
					for(int m = 0;m < NUM_MUTATION;++m){
						baseline[0](c, r, m, w) = mutation_freq(c, r, m, w);
					}
				}
			}
		}

		baseline_scores = predict_mutations(baseline, var, test, db, diff, regions,
			opt.motif_len);

		ave_baseline_score = average_first( baseline_scores.begin(), baseline_scores.end() );

		baseline_scores.clear();

		if(IS_ROOT){
			cerr << "\tTraining distribution baseline model score = " << ave_baseline_score << endl;
		}
		
		#ifdef GRAPHIC_FOR_1663
		// DEBUG -- for the LANL 1663 graphic illustrating model output at each base position
		// Use the "reference" genome as the parent

		mutation_freq.normalize();
		
		const Sequence &ref_seq = db["reference"];

		for(Sequence::const_iterator r = ref_seq.begin();r != ref_seq.end();++r){
			cout << '\t' << r->ascii_base();
		}

		cout << endl;

		for(unsigned int m = 0;m < var.num_mutation;++m){

			switch(m){
				case MUTATION_SUB_A:
					cout << ">A";
					break;
				case MUTATION_SUB_T:
					cout << ">T";
					break;
				case MUTATION_SUB_G:
					cout << ">G";
					break;
				case MUTATION_SUB_C:
					cout << ">C";
					break;
				case MUTATION_INSERTION:
					cout << "ins";
					break;
				case MUTATION_DELETION:
					cout << "del";
					break;
				default:
					throw __FILE__ ":main: Unknown mutation type";
			};

			for(unsigned int i = 0;i < opt.motif_len/2;++i){
				cout << '\t';
			}

			ForEachSenseWord(ref_seq, opt.motif_len)
		
				float score = -1.0;
				bool is_blank = true;

				// Compute scores for full length words
				if(ValidWord){

					// The coordinates of the center of the word (in the parent frame of reference)
					const unsigned int center_loc = Loc3 - opt.motif_len/2;

					// Extract the center base by shifting two bits per base for half of the motif
					// length
					const unsigned char center_base = ( SenseWord >> 2*(opt.motif_len/2) ) & 3;

					// Evaluate the predicted mutation scores for the base in the center of this word
					const unsigned int region_index = 0; // Force region to zero

					// Lookup the precomputed score
					score = mutation_freq(ref_seq[center_loc].codon(), region_index, 
						m, SenseWord);

					is_blank = false;

					if( (m == MUTATION_SUB_A) && (center_base == BASE_A) ){
						is_blank = true;
					}

					if( (m == MUTATION_SUB_T) && (center_base == BASE_T) ){
						is_blank = true;
					}

					if( (m == MUTATION_SUB_G) && (center_base == BASE_G) ){
						is_blank = true;
					}

					if( (m == MUTATION_SUB_C) && (center_base == BASE_C) ){
						is_blank = true;
					}

					cout << '\t';
				
					if(!is_blank){
						cout << score;
					}
				}

			EndWord

			cout << endl;
		}

		throw __FILE__ ":main: Output graphic for LANL 1663 complete!";
		#endif // GRAPHIC_FOR_1663

		// Evaluate any models that the user has provided
		for(deque<string>::const_iterator i = opt.eval_filename.begin();i != opt.eval_filename.end();++i){

			cerr << "Evaluating model parameters from " << *i << endl;

			if(IS_ROOT){
				baseline[0] = load_param(*i, var, num_word);
			}

			broadcast(baseline, mpi_rank, 0);

			baseline_scores = predict_mutations(baseline, var, test, db, diff, regions,
				opt.motif_len);

			ave_baseline_score = average_first( baseline_scores.begin(), baseline_scores.end() );

			baseline_scores.clear();

			if(IS_ROOT){
				cerr << "\tTraining distribution score = " << ave_baseline_score << endl;
			}
		}

		// Motif (kmer) are indexed as:
		//		motif[ensemble index][codon_loc][region][redundant motif][mutation type]
		// The first index, ensemble member, is stored as a separate deque
		deque< Tensor<float> > epoch_param;

		for(size_t epoch = 0;epoch < opt.num_epoch;++epoch){

			// Randomize the training set
			randomize(train.begin(), train.end(), rand_gen);

			// Make sure that all MPI ranks agree on the ordering of train by using the
			// ordering determined by rank 0
			broadcast(train, mpi_rank, 0);

			for(size_t batch = 0;batch < opt.num_batch;++batch){

				if(IS_ROOT){
					cerr << "Processing batch " << (batch + 1) << "/" << opt.num_batch << "; epoch " 
						<< (epoch + 1) << "/" << opt.num_epoch << "; fitting "
						<< opt.motif_len << " bp motifs" << endl;
				}

				deque< pair<string, string> > local_train;

				for(size_t j = 0;j < num_train;++j){

					if(j%opt.num_batch == batch){
						local_train.push_back(train[j]);
					}
				}

				if(IS_ROOT){
					cerr << "\t|train| = " << local_train.size() << endl;
				}

				epoch_param.push_back( find_motifs(var, local_train, db, diff, regions, opt, rand_gen) );
			}
		}

		if( epoch_param.empty() ){
			throw __FILE__ ":main: |epoch_param| = 0";
		}

		deque< pair<float, string> > test_scores = predict_mutations(epoch_param, var, test, db, diff, regions,
			opt.motif_len);

		const float test_score = average_first( test_scores.begin(), test_scores.end() );

		if(IS_ROOT){
			
			cerr << "For k = " << opt.motif_len << "; final test set score = " 
				<< test_score << endl;

			const size_t ensemble_size = epoch_param.size();

			cerr << "k-mer ensemble size = " << ensemble_size << endl;

			for(size_t sample = 0;sample < ensemble_size;++sample){

				cerr << "############## Ensemble sample[" << sample << "] ##############" << endl;

				for(unsigned int m = 0;m < var.num_mutation;++m){

					for(unsigned int c = 0;c < var.num_codon_loc;++c){

						for(unsigned int r = 0;r < var.num_region;++r){
	
							cerr << "k-mer motif for mutation " << mutation_name[m] 
								<< "; codon loc " << c 
								<< "; region " << r
								<< endl;

							for(Word w = 0;w < num_word;++w){
								cerr << word_to_string(w, opt.motif_len) << '\t' 
									<< epoch_param[sample](c, r, m, w) << endl;
							}
						}
					}
				}

				cerr << endl;
			}

			if(false){

				// Print the per-sequence test scores
				sort( test_scores.begin(), test_scores.end() );

				cerr << "Per-sequence test AUROC values:" << endl;

				// Use descending order
				for(deque< pair<float, string> >::const_reverse_iterator i = test_scores.rbegin();i != test_scores.rend();++i){
					cerr << i->second << '\t' << i->first << endl;
				}
			}
		}

		gsl_rng_free(rand_gen);

		profile = MPI_Wtime() - profile;

		if(IS_ROOT){
			cerr << "Model fitting and evaluation complete in " << profile << " sec" << endl;
		}

		MPI_Finalize();
	}
	catch(const char *error){

		cerr << "[" << mpi_rank << "] Caught the error: " << error << endl;

		MPI_Finalize();

		return EXIT_FAILURE;
	}
	catch(const exception &error)
	{
		cerr << "[" << mpi_rank << "] Caught the error: " << error.what() << endl;

		MPI_Finalize();

		return EXIT_FAILURE;
	}
	catch(...){

		cerr << "[" << mpi_rank << "] Caught an unhandled error" << endl;

		MPI_Finalize();

		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

ostream& operator<<(ostream& m_out, const Motif& m_motif)
{
	const unsigned int len = m_motif.size();

	for(int base = 0;base < NUM_BASE;++base){
		
		if(base != 0){
			m_out << '\n';
		}

		m_out << base_name[base];
		
		for(unsigned int i = 0;i < len;++i){
			m_out << '\t' << m_motif(base, i);
		}
	}

	return m_out;
}

void load_genealogy(MULTIMAP<string /*parent*/, string /*child*/> &m_family, MAP<string /*child*/, GenomeDifference> &m_diff,
	const string &m_filename)
{
	ifstream fin( m_filename.c_str() );

	if(!fin){
		throw __FILE__ ":load_genealogy: Unable to open input file";
	}

	string line;

	size_t skipped = 0;

	// Parse the three-column, tab-delimited text file: [child accession][parent accession][difference string]
	while( getline(fin, line) ){

		// Skip text after a comment
		string::size_type loc = line.find('#');

		if(loc != string::npos){

			// Trim the comment text
			line = line.substr(0, loc);
		}

		// Skip blank lines
		if( line.empty() ){
			continue;
		}

		stringstream buffer(line);

		string parent;
		string child;
		string difference;

		if( !(buffer >> child) ){
			throw __FILE__ ":load_genealogy: Unable to parse child accession";
		}

		if( !(buffer >> parent) ){
			throw __FILE__ ":load_genealogy: Unable to parse parent accession";
		}

		// Don't test for successful reading of the difference string, as this string is allowed to be empty
		// (i.e., the parent and child are identical).
		buffer >> difference;

		// Don't store parent-child pairs unless there is a difference between the parent and child
		if( !difference.empty() ){

			m_family.insert( make_pair(parent, child) );

			// Make sure that each child is descended from one, and only one, parent
			if( m_diff.find(child) != m_diff.end() ){
				throw __FILE__ ":load_genealogy: Multiple parents for the same child";
			}

			m_diff[child] = GenomeDifference(difference);
		}
		else{
			++skipped;
		}
	}

	if(IS_ROOT){
		cerr << "Skipped " << skipped << " identical parent-child pairs" << endl;
	}
}

// Compute a motif for each NUM_MUTATION possible mutation outcome:
// 		Central base -> A
// 		Central base -> T
// 		Central base -> G
// 		Central base -> C
// 		Central base -> single base deletion
// 		Central base -> One or more base insertion
//
//	All of these k-mer motifs are computed at the same time, so that the highest-scoring
// 	k-mer motif corresponds to actual mutations observed in the training data.
Tensor<float> find_motifs(const TensorDim &m_var, const deque< pair<string, string> > &m_train, 
	const MAP<string, Sequence> &m_db,
	const MAP<string, GenomeDifference> &m_diff, const vector< pair<unsigned int, unsigned char> > &m_regions,
	const MotifOptions &m_opt, const gsl_rng *m_rand_ptr)
{
	// Make sure that the user has specified an *odd* motif length (so we can determine the central base)
	if(m_opt.motif_len%2 != 1){
		throw __FILE__ ":find_motifs: Please specify an odd motif length";
	}

	const unsigned int num_word = 1 << 2*m_opt.motif_len;
	
	// The objective function that will be optimized to determine the k-mer motifs is the sum of accuracy values for
	// each element in the training set.
	// 		- For each parent-child pair, there is a set of one or more *observed* mutations (consisting of a mutation type
	//		  occurring at a specific location).
	//		- A k-mer motif maps to a quantitative score associated with each mutation type at each location.
	//		- The objective function should reward k-mer motifs that map to a higher score for observered mutations than for other
	//		  possible mutations.
	//		- Computing the scores for all possible mutations yeilds (approximately) NUM_MUTATION x (genome length) scalar values.
	//		- Each of these scalar values is associated with either an observed mutation, "O", or a potential mutation that is *not*
	//		  observed for this pair, "p".
	//		- For each parent-child pair, sort the motif scores (and associate score labels) in ascending order and compute the 
	//		  maximum possible accuracy (allowing for a varrying "threshold" score):
	//				pppppOppppppp OOOOOO
	//							 ^
	//						Best accuracy threshold (Tp = 6, Tn = 12, Fp = 0, Fn = 1) --> accuracy = (6 + 12)/(6 + 12 + 0 + 1) = 0.947

	const unsigned int num_train = m_train.size();
	
	// Define a default region for the case when the user does not specify *any* regions
	const unsigned int default_region = 0;

	// For a motif of length m_opt.motif_len, there are 4^m_opt.motif_len number of possible motifs
	
	// Step 1: Find the set of observed (and potential) mutations for each pair of genome sequences in the training set.
	vector< deque<WordIndex> > observed(num_train);

	if(IS_ROOT){
		cerr << "\tAssembling training fragments: ";
	}

	string info;

	const int update_every = max(1U, num_train/10);

	for(unsigned int i = 0;i < num_train;++i){

		if( IS_ROOT && (i%update_every == 0) ){

			stringstream ssin;

			ssin << (100.0*i)/num_train << '%';

			for(string::const_iterator j = info.begin();j != info.end();++j){
				cerr << '\b';
			}

			for(string::const_iterator j = info.begin();j != info.end();++j){
				cerr << ' ';
			}

			for(string::const_iterator j = info.begin();j != info.end();++j){
				cerr << '\b';
			}

			info = ssin.str();

			cerr << info;
		}

		MAP<string /*accession*/, Sequence>::const_iterator parent_iter = m_db.find(m_train[i].first);

		if( parent_iter == m_db.end() ){
			throw __FILE__ ":find_motifs: Unable to lookup parent sequence";
		}

		MAP<string, GenomeDifference>::const_iterator child_iter = m_diff.find(m_train[i].second);

		if( child_iter == m_diff.end() ){
			throw __FILE__ ":find_motifs: Unable to lookup child GenomeDifference";
		}

		deque< pair<unsigned int, char> >::const_iterator sub_iter = child_iter->second.get_substitution().begin();
		deque< pair<unsigned int, string> >::const_iterator ins_iter = child_iter->second.get_insertion().begin();
		deque<unsigned int>::const_iterator del_iter = child_iter->second.get_deletion().begin();

		const deque< pair<unsigned int, char> >::const_iterator sub_end = child_iter->second.get_substitution().end();
		const deque< pair<unsigned int, string> >::const_iterator ins_end = child_iter->second.get_insertion().end();
		const deque<unsigned int>::const_iterator del_end = child_iter->second.get_deletion().end();

		vector< pair<unsigned int, unsigned char> >::const_iterator region_iter = m_regions.begin();
		vector< pair<unsigned int, unsigned char> >::const_iterator next_region_iter = region_iter;

		// Regions are optional
		if( next_region_iter != m_regions.end() ){
			++next_region_iter;
		}

		// A buffer to count the number of occurances of each word in each region, codon position and mutation state
		Tensor<unsigned int> buffer(m_var.num_codon_loc, 
			m_var.num_region, 
			m_var.num_mutation + 1, // Add 1 to include the "no mutation" state
			num_word);

		ForEachSenseWord(parent_iter->second, m_opt.motif_len)
		
			// Only store full length words
			if(ValidWord){

				// The coordinates of the center of the word (in the parent frame of reference)
				const unsigned int center_loc = Loc3 - m_opt.motif_len/2;

				// Advance the region if needed
				if( ( next_region_iter != m_regions.end() ) && (next_region_iter->first <= center_loc) ){

					region_iter = next_region_iter;
					++next_region_iter;
				}

				unsigned char mutation_state = NO_MUTATION;

				// Advance the iterators (if needed)
				while( (sub_iter != sub_end) && (sub_iter->first < center_loc) ){
					++sub_iter;
				}

				while( (ins_iter != ins_end) && (ins_iter->first < center_loc) ){
					++ins_iter;
				}

				while( (del_iter != del_end) && (*del_iter < center_loc) ){
					++del_iter;
				}

				if( (sub_iter != sub_end) && (sub_iter->first == center_loc) ){
					
					// This parent fragment maps to a substitution in the child
					if(mutation_state != NO_MUTATION){
						throw __FILE__ ":find_motifs: Found overlapping mutation states";
					}

					switch(sub_iter->second){
						case 'A': case 'a':
							mutation_state = MUTATION_SUB_A;
							break;
						case 'T': case 't':
							mutation_state = MUTATION_SUB_T;
							break;
						case 'G': case 'g':
							mutation_state = MUTATION_SUB_G;
							break;
						case 'C': case 'c':
							mutation_state = MUTATION_SUB_C;
							break;
					};

					++sub_iter;
				}

				#ifdef MUTATION_INSERTION
				if( (ins_iter != ins_end) && (ins_iter->first == center_loc) ){
					
					// This parent fragment maps to an insertion in the child
					if(mutation_state != NO_MUTATION){
						throw __FILE__ ":find_motifs: Found overlapping mutation states";
					}

					mutation_state = MUTATION_INSERTION;

					++ins_iter;
				}
				#endif // MUTATION_INSERTION

				#ifdef MUTATION_DELETION
				if( (del_iter != del_end) && (*del_iter == center_loc) ){
					
					// This parent fragment maps to a deletion in the child
					if(mutation_state != NO_MUTATION){
						throw __FILE__ ":find_motifs: Found overlapping mutation states";
					}

					mutation_state = MUTATION_DELETION;

					++del_iter;
				}
				#endif // MUTATION_DELETION

				// Count the number of times we observe a given word associated with
				// a given mutation state
				const unsigned int region_index = ( region_iter != m_regions.end() ) ? region_iter->second : default_region;

				// DEBUG
				//cerr << "region_index for base " << center_loc << " is " << region_index << endl;
				
				buffer(parent_iter->second[center_loc].codon(), region_index, mutation_state, SenseWord) += 1;
			}

		EndWord

		deque<WordIndex> &obs = observed[i];

		for(unsigned int c = 0;c < m_var.num_codon_loc;++c){
			
			for(unsigned int r = 0;r < m_var.num_region;++r){

				// Include the "no mutation" state
				for(unsigned int m = 0;m < m_var.num_mutation + 1;++m){
				
					for(unsigned int w = 0;w < num_word;++w){

						unsigned int count = buffer(c, r, m, w);

						if(count > 0){

							obs.push_back( 
								WordIndex(
									count, // Weight
									w, // numeric value of kmer word
									r, // region of the genome
									c, // codon location of central base
									m, // mutation state
									( w>> 2*(m_opt.motif_len/2) ) & 3) // center base
								);
						}
					}
				}
			}
		}

		// To avoid clobbering the CPU cache, sort the "observed" words so that mutations that are 
		// actually observed are at the begining of the deque.
		sort(obs.begin(), obs.end(), sort_by_mutation());
	}

	if(IS_ROOT){
		cerr << " complete!" << endl;
	}

	// Count the number of mutation states, not including the "no mutation" state
	vector<size_t> count(m_var.num_mutation);

	for(vector< deque<WordIndex> >::const_iterator i = observed.begin();i != observed.end();++i){

		for(deque<WordIndex>::const_iterator j = i->begin();j != i->end();++j){

			if(j->state < m_var.num_mutation){
				count[j->state] += j->weight;
			}
		}
	}

	if(IS_ROOT){

		cerr << "\tMutation inventory: ";
		cerr << count[MUTATION_SUB_A] << " -> A";
		cerr << ", " << count[MUTATION_SUB_T] << " -> T";
		cerr << ", " << count[MUTATION_SUB_G] << " -> G";
		cerr << ", " << count[MUTATION_SUB_C] << " -> C";

		#ifdef MUTATION_DELETION
			cerr << ", " << count[MUTATION_DELETION] << " del";
		#endif // MUTATION_DELETION

		#ifdef MUTATION_INSERTION
			cerr << ", " << count[MUTATION_INSERTION] << " ins";
		#endif // MUTATION_INSERTION

		cerr << endl;
	}

	info = "";

	if(IS_ROOT){
		cerr << "\tOptimizing:" << endl;
	}

	// Greedy optimization method #1
	Tensor<float> best_motif(m_var.num_codon_loc, m_var.num_region, m_var.num_mutation, num_word);

	// Initialize all motif values to zero
	best_motif = 0.0;

	float best_score = score(observed, best_motif, m_opt.motif_len);

	// Start in the center of each motif, test each base, and radiate out...

	// Keep iterating until the score stops improving AND we've reduced the step size
	enum {
		IMPROVED, 
		NOT_IMPROVED, 
		UNABLE_TO_IMPROVE
	};

	int search_status = IMPROVED;

	int num_step = 4;

	for(unsigned int iteration = 0;search_status != UNABLE_TO_IMPROVE;++iteration){

		double profile = MPI_Wtime();

		Tensor<float> local_best_motif(best_motif);
		
		++search_status;

		const float delta_step = 1.0/num_step;

		for(unsigned int curr_word = 0;curr_word < num_word;++curr_word){
			
			const Word central_base = ( curr_word >> 2*(m_opt.motif_len/2) ) & 3;

			for(unsigned int r = 0;r < m_var.num_region;++r){

				for(unsigned int c = 0;c < m_var.num_codon_loc;++c){

					for(unsigned int m = 0;m < m_var.num_mutation;++m){

						if(m == central_base){
							// Self-substitutions are not allowed, so skip the searching for central
							// positions that have a base substituting to itself
							continue;
						}

						// Make sure all ranks agree on the same best motif and score
						//broadcast(best_score, mpi_rank, 0);

						// Save a copy of the specific word value so we can restore the state of best_motif
						const float original_x = best_motif(c, r, m, curr_word);

						// Test k-mer element values in descending order
						for(int step = num_step;step >= 0;--step){
							
							// x -> [-1, 1]
							const float x = 2.0*step*delta_step - 1.0;

							best_motif(c, r, m, curr_word) = x;

							const float s = score(observed, best_motif, m_opt.motif_len);

							// DEBUG
							//if(IS_ROOT){
							//	cerr << "+" << word_to_string(curr_word, m_opt.motif_len) 
							//		<< '\t' << mutation_name[m] << '\t' 
							//		<< x << '\t' << s << endl;
							//}

							if(s > best_score){ // The first (highest) element value
								
								best_score = s;
								local_best_motif = best_motif;
								
								search_status = IMPROVED;

								if(IS_ROOT){
									cerr << word_to_string(curr_word, m_opt.motif_len) 
										<< '\t' << mutation_name[m] << '\t' 
										<< x << '\t' << s;
									
									if(m_var.num_codon_loc > 1){
										cerr << "\t(codon " << c << ")";
									}

									if(m_var.num_region > 1){
										cerr << "\t(region " << r << ")";
									}

									cerr << endl;
								}
							}

							//cout << base_name[b] << '\t' << j << '\t' << mutation_name[m] << '\t' 
							//	<< x*delta << '\t' << s << " (" << best_score << ")" << endl;
						}

						// Restore the state of the best_motif
						best_motif(c, r, m, curr_word) = original_x;

						//cout << endl;
					}
				}
			}
		}

		if(IS_ROOT){

			profile = MPI_Wtime() - profile;

			cerr << iteration << '\t' << best_score << " in " << profile << " sec" << endl;
		}

		if(search_status == IMPROVED){
			best_motif = local_best_motif;
		}

		if(search_status == NOT_IMPROVED){
			
			//num_step *= 4;
			num_step *= 2;
			
			if(IS_ROOT){
				cerr << "Increased number of steps to " << num_step << endl;
			}
		}

		// Make sure all ranks agree on the same best motif and score
		//broadcast(best_motif, mpi_rank, 0);
	}

	if(IS_ROOT){
		cerr << endl;
	}

	return best_motif;
}

// m_mut -> |number of parents| x |parent word|
// m_buffer -> tensor of cached score values, with index order:
//		codon_loc
//		region
//		mutatation
//		word
float score(const vector< deque<WordIndex> > &m_mut, const Tensor<float> &m_buffer, const unsigned int &m_kmer_len)
{
	float ret = 0.0;

	// Count the number of parent-child pairs with one or more observed mutations
	unsigned int num_valid_parent_child = 0;

	const int num_mut = m_mut.size();

	// Iterate over training pairs. Each MPI rank processes a unique subset of parent-child pairs
	#pragma omp parallel for \
			reduction(+:ret, num_valid_parent_child)
	for(int i = mpi_rank;i < num_mut;i += mpi_numtasks){

		const deque<WordIndex>& words = m_mut[i];

		// Due to the unbalanced nature of the problem, we expect very few observed mutations!
		// Update -- the imbalanced data problem is more subtle than orginally thought. There appears
		// to be a power-law distribution of mutation counts (i.e., a large number of parent-child pairs
		// that only differ by a small number of mutations and a small number of parent-child pairs
		// that differ by a large number of mutations).
		deque<float> obs; // Observed mutations -- will be either true positive or false negative

		// Lookup the observed mutation scores for the current genome.
		// The words deque has been sorted, so that observed mutations are first. When we
		// get to the first NO_MUTATION state, then we have exhausted the observed mutations.
		for(deque<WordIndex>::const_iterator w = words.begin();( w != words.end() ) && (w->state != NO_MUTATION);++w){

			#ifdef EBUG
			if(w->state > NUM_MUTATION){
				throw __FILE__ ":score: Invalid state detected!";
			}
			#endif // EBUG

			// The weight is the count of the number of times a given motif kmer and
			// associated mutation state is found in a genome. For observed mutations,
			// we expect the weight to be small (most likely < 5?), so just replicate
			// the score by pushing back multiple times.
			const float s = m_buffer(w->codon_loc, w->region, w->state, w->index);

			for(unsigned int j = 0;j < w->weight;++j){
				obs.push_back(s);
			}
		}

		// Skip training pairs with no observed mutations
		if( obs.empty() ){
			continue;
		}

		// Observation: About %50 of the time is spent evaluating the Mann-Whitney score for
		// parent-child pairs that have more than five mutations. How can we speed up this
		// calculation?
		//	- Sorting the observed scores (in obs) in order to break out early does not
		//	  appear to speed up execution.

		++num_valid_parent_child;

		// Count the number of potential (i.e. unobserved) mutation scores that are less than,
		// or tied with, each of the observed scores
		unsigned int num_tie = 0U;
		unsigned int num_less = 0U;
		unsigned int total_weight = 0U;

		//float total_rank = 0.0;

		// Take advantage of the imbalance between observed and unobserved mutations. Stream the 
		// unobserved mutation scores and count the number which are less than each of the observed
		// mutation scores. This avoids a very expensive sorting operation.
		for(deque<WordIndex>::const_iterator w = words.begin();w != words.end();++w){

			// Extract the center base by shifting two bits per base for half of the motif
			// length
			//const unsigned char center_base = ( w->index >> 2*(m_kmer_len/2) ) & 3;

			// For each word, compute the maximum predicted mutation score over all mutation predictions
			//float max_pot_score = std::numeric_limits<float>::lowest();

			for(unsigned char m = 0;m < NUM_MUTATION;++m){

				// 1) Do not include the score of an observed mutations
				// 2) Forbid substitutions of a base for itself
				if( (w->center_base != m) && (w->state != m) ){

					const float score = m_buffer(w->codon_loc, w->region, m, w->index);

					// Note that w->weight is the number of times a particular word and an associated mutation
					// state is observed in a given genome sequence.
					total_weight += w->weight;

					for(deque<float>::const_iterator k = obs.begin();k != obs.end();++k){
					
						// Count the number of unobserved mutations that are *less* than this observed mutation score.
						// Follow the Mann-Whitney U test rule of counting ties as 0.5
						// (see https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test).
						const float delta = score - *k;

						num_less += (delta < -EPS) ? w->weight : 0U;
						num_tie += (fabsf(delta) < EPS) ? w->weight : 0U;
					}
				}
			}
		}

		const float total_rank = num_tie*0.5f + num_less;

		// Accumulate the average *fractional* rank of the observed mutations (i.e. positive examples)
		ret += !words.empty() ? (float)(total_rank)/( total_weight*obs.size() ) : 0; 
	}

	// Accumulate the average fraction rank across all MPI ranks
	if(MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS){
		throw __FILE__ ":score: Error in MPI_Allreduce(ret)";
	}

	// Accumulate the number of valid parent-child pairs across all MPI ranks
	if(MPI_Allreduce(MPI_IN_PLACE, &num_valid_parent_child, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS){
		throw __FILE__ ":score: Error in MPI_Allreduce (num_valid_parent_child)";
	}

	// Return the average (over training pairs) of the average fractional rank
	if(num_valid_parent_child > 0){
		ret /= num_valid_parent_child;
	}

	return ret;
}

deque< pair<float /*child accession*/, string /*child auroc*/> > 
	predict_mutations(const deque< Tensor<float> > &m_param, const TensorDim &m_var,
	const deque< pair<string, string> > &m_test, const MAP<string, Sequence> &m_db, 
	const MAP<string, GenomeDifference> &m_diff, const vector< pair<unsigned int, unsigned char> > &m_regions,
	const unsigned int &m_motif_len)
{
	const size_t num_ensemble = m_param.size();

	if(num_ensemble == 0){
		throw __FILE__ ":predict_mutations: |m_param| = 0";
	}

	if(IS_ROOT){
		cerr << "\tPredicting mutations for " << m_test.size() << " parent-child pairs with " << num_ensemble
			<< " ensemble models" << endl;
	}
	
	// Use the actual motif_len, not opt.motif_len
	if( m_param.front().empty() ){
		throw __FILE__ ":predict_mutations: Empty tensor";
	}

	//const unsigned int num_word = 1 << 2*motif_len;

	//////////////////////////////////////////////////////////////////
	// We have an ensemble of models that we need to evaluate against
	// each word and compute the average score. Enumerate this average
	// score of all possible input sequences
	
	// Make sure that all ensemble models have the *same* motif length (this is required if we want to
	// precompute the motif sequence to motif score lookup table)
	for(deque< Tensor<float> >::const_iterator param_iter = m_param.begin();param_iter != m_param.end();++param_iter){
		
		if( param_iter->empty() ){
			throw __FILE__ ":predict_mutations: |*param_iter| = 0 ";
		}

		// Make sure we have consistent tensor indexing
		if( m_var.num_codon_loc != param_iter->size(0) ){
			throw __FILE__ ":predict_mutations: Unexpected size of num_codon_loc";
		}

		if( m_var.num_region != param_iter->size(1) ){
			throw __FILE__ ":predict_mutations: Unexpected size of num_region";
		}

		if( m_var.num_mutation != param_iter->size(2) ){
			throw __FILE__ ":predict_mutations: Unexpected size of num_mutation";
		}
	}

	deque< pair<float, string> > auroc;

	const int num_test = m_test.size();

	// Each MPI rank evalutes the scores for a unique subset of parent-child pairs
	#pragma omp parallel
	{
		deque< pair<float, string> > local_auroc;

		#pragma omp for
		for(int i = mpi_rank;i < num_test;i += mpi_numtasks){
			
			// DEBUG
			//if(mpi_rank == 1){
			//	cout << i << '\t' << m_test[i].first << '\t' << m_test[i].second << endl;
			//}

			MAP<string, Sequence>::const_iterator parent_iter = m_db.find(m_test[i].first);

			if( parent_iter == m_db.end() ){
				throw __FILE__ ":predict_mutations: Unable to lookup parent sequence";
			}

			MAP<string, GenomeDifference>::const_iterator child_iter = m_diff.find(m_test[i].second);

			if( child_iter == m_diff.end() ){
				throw __FILE__ ":predict_mutations: Unable to lookup child GenomeDifference";
			}

			deque< pair<unsigned int, char> >::const_iterator sub_iter = child_iter->second.get_substitution().begin();
			deque< pair<unsigned int, string> >::const_iterator ins_iter = child_iter->second.get_insertion().begin();
			deque<unsigned int>::const_iterator del_iter = child_iter->second.get_deletion().begin();

			const deque< pair<unsigned int, char> >::const_iterator sub_end = child_iter->second.get_substitution().end();
			const deque< pair<unsigned int, string> >::const_iterator ins_end = child_iter->second.get_insertion().end();
			const deque<unsigned int>::const_iterator del_end = child_iter->second.get_deletion().end();

			vector< pair<unsigned int, unsigned char> >::const_iterator region_iter = m_regions.begin();
			vector< pair<unsigned int, unsigned char> >::const_iterator next_region_iter = region_iter;

			// Regions are optional
			if( next_region_iter != m_regions.end() ){
				++next_region_iter;
			}

			// The scores for the observed mutations in the child genome
			MULTIMAP<int /*mutation*/, float /*score*/> obs_score;

			// The scores for potential, but unobserved, mutations in the child genome
			//deque<float> pot_score;
			deque< pair<int /*mutation*/, float> > pot_score;

			ForEachSenseWord(parent_iter->second, m_motif_len)
			
				// Compute scores for full length words
				if(ValidWord){

					// The coordinates of the center of the word (in the parent frame of reference)
					const unsigned int center_loc = Loc3 - m_motif_len/2;

					// Extract the center base by shifting two bits per base for half of the motif
					// length
					const unsigned char center_base = ( SenseWord >> 2*(m_motif_len/2) ) & 3;

					// Advance the region if needed
					if( ( next_region_iter != m_regions.end() ) && (next_region_iter->first <= center_loc) ){

						region_iter = next_region_iter;
						++next_region_iter;
					}

					unsigned char mutation_state = NO_MUTATION;

					// Advance the child iterators (if needed)
					while( (sub_iter != sub_end) && (sub_iter->first < center_loc) ){
						++sub_iter;
					}

					while( (ins_iter != ins_end) && (ins_iter->first < center_loc) ){
						++ins_iter;
					}

					while( (del_iter != del_end) && (*del_iter < center_loc) ){
						++del_iter;
					}

					if( (sub_iter != sub_end) && (sub_iter->first == center_loc) ){
						
						// This parent fragment maps to a substitution in the child
						if(mutation_state != NO_MUTATION){
							throw __FILE__ ":predict_mutations: Found overlapping mutation states";
						}

						switch(sub_iter->second){
							case 'A': case 'a':
								mutation_state = MUTATION_SUB_A;
								break;
							case 'T': case 't':
								mutation_state = MUTATION_SUB_T;
								break;
							case 'G': case 'g':
								mutation_state = MUTATION_SUB_G;
								break;
							case 'C': case 'c':
								mutation_state = MUTATION_SUB_C;
								break;
						};

						++sub_iter;
					}

					#ifdef MUTATION_INSERTION
					if( (ins_iter != ins_end) && (ins_iter->first == center_loc) ){
						
						// This parent fragment maps to an insertion in the child
						if(mutation_state != NO_MUTATION){
							throw __FILE__ ":predict_mutations: Found overlapping mutation states";
						}

						mutation_state = MUTATION_INSERTION;

						++ins_iter;
					}
					#endif //MUTATION_INSERTION

					#ifdef MUTATION_DELETION
					if( (del_iter != del_end) && (*del_iter == center_loc) ){
						
						// This parent fragment maps to a deletion in the child
						if(mutation_state != NO_MUTATION){
							throw __FILE__ ":predict_mutations: Found overlapping mutation states";
						}

						mutation_state = MUTATION_DELETION;

						++del_iter;
					}
					#endif // MUTATION_DELETION

					// Evaluate the predicted mutation scores for the base in the center of this word
					const unsigned int region_index = ( region_iter != m_regions.end() ) ?
						region_iter->second : 0;

					for(unsigned int m = 0;m < m_var.num_mutation;++m){

						// Lookup the precomputed score
						const float local = m_param.front()(parent_iter->second[center_loc].codon(), region_index, 
							m, SenseWord);

						if(m == mutation_state){
							obs_score.insert( make_pair(mutation_state, local) );
						}
						else{
							
							// Record the predicted mutation score. Bases are *not*
							// allowed to substitute for themselves!
							if(m != center_base ){
								pot_score.push_back( make_pair(m, local) );
							}
						}
					}
				}

			EndWord

			// We expect that |obs_score| + |pot_score| ~ number of valid bases in a genome
			if( obs_score.empty() ){

				// Skip child genomes that don't have any observable mutational difference from the parent
				// (could be due to 'N' or other degenerate bases in the parent)
				continue;
			}

			if( pot_score.empty() ){
				throw __FILE__ ":predict_mutations: Did not find any potential (but unobserved) mutations!";
			}

			// Compute the average, fractional rank of observed mutation scores for this child genome
			unsigned int total_num_tie = 0U;
			unsigned int total_num_less = 0U;

			for(deque< pair<int, float> >::const_iterator p = pot_score.begin();p != pot_score.end();++p){

				for(MULTIMAP<int, float>::const_iterator k = obs_score.begin();k != obs_score.end();++k){

					// Follow the Mann-Whitney U test rule of counting ties as 0.5
					// (see https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test).
					//total_rank += (pot_score[j] < *k);
					const float delta = p->second - k->second;

					total_num_tie += (fabs(delta) < EPS) ? 1U : 0U;
					total_num_less += (delta < -EPS)  ? 1U : 0U;

					// DEBUG
					//if(delta < -EPS){

					//	cerr << delta << '\t' << k->second << " (" << mutation_name[k->first] << "); " 
					//		<< p->second << " (" << mutation_name[p->first] << ")" << endl;
					//}
				}
			}

			const float total_rank = 0.5*total_num_tie + total_num_less;

			//auroc += total_rank/( pot_score.size() * obs_score.size() );
			local_auroc.push_back( make_pair(
				total_rank/( pot_score.size() * obs_score.size() ), /*per-child AUROC*/
				m_test[i].second /*child accession*/
			) );
		}

		// Merge the thread-local AUROC scores
		#pragma omp critical
		while( !local_auroc.empty() ){

			auroc.push_back( local_auroc.back() );
			local_auroc.pop_back();
		}
	}

	// DEBUG
	//cerr << "[" << mpi_rank << "] |auroc| = " << auroc.size() << "; |test| = " << m_test.size() << endl;

	// Gather the per-child AUROC values from all ranks
	deque< pair<float, string> > global_auroc;

	for(int i = 0;i < mpi_numtasks;++i){

		deque< pair<float, string> > tmp(auroc);

		broadcast(tmp, mpi_rank, i);

		while( !tmp.empty() ){

			global_auroc.push_back( tmp.back() );
			tmp.pop_back();
		}
	}

	// DEBUG
	//if(mpi_rank == 0){
		//for(deque< pair<float, string> >::const_iterator i = global_auroc.begin();i != global_auroc.end();++i){
		//	cout << i->second << '\t' << i->first << endl;
		//}

		//for(deque< pair<string, string> >::const_iterator i = m_test.begin();i != m_test.end();++i){
		//	cout << i->first << '\t' << i->second << endl;
		//}
	//}

	//throw "DEBUG";

	return global_auroc;
}

// Approximate identification of potential secondary structure elements by searching for
// complementary sequences
vector<bool> mask_secondary_structure(const Sequence &m_seq, const unsigned int &m_kmer_len,
	const unsigned int &m_window)
{

	vector<bool> mask(m_seq.size(), false);

	MULTIMAP<Word, int> seq_loc;

	ForEachSenseWord(m_seq, m_kmer_len)

		if(ValidWord){
			seq_loc.insert( make_pair(SenseWord, Loc5) );
		}

	EndWord

	//sort( seq_set.begin(), seq_set.end() );
	//seq_set.erase( unique( seq_set.begin(), seq_set.end() ), seq_set.end() );

	//unsigned int count = 0;

	ForEachAntisenseWord(m_seq, m_kmer_len)

		if(ValidWord){
			
			typedef MULTIMAP<Word, int>::const_iterator I;

			const pair<I, I> range = seq_loc.equal_range(AntisenseWord);

			for(I i = range.first;i != range.second;++i){
				
				if( abs( i->second - int(Loc5) ) < m_window){
					
					//++count;

					for(unsigned int j = Loc5;j <= Loc3;++j){
						mask[j] = true;
					}
				}
			}
		}

	EndWord

	// How many bp are masked?
	//unsigned int num_base = 0;

	//for(vector<bool>::const_iterator i = mask.begin();i != mask.end();++i){
	//	num_base += *i;
	//}

	//cerr << "Found " << count << " complementary sequences masking " << num_base << " bp/" << mask.size() << " bp" << endl;

	return mask;
}

// Recursively compute the distance to the root node
size_t distance_to_root(const string &m_node, const MAP<string, string> &m_child_to_parent)
{
	MAP<string, string>::const_iterator i = m_child_to_parent.find(m_node);

	if( i == m_child_to_parent.end() ){
		return 0;
	}

	return 1 + distance_to_root(i->second, m_child_to_parent);
}

unsigned char get_region(const vector< pair<unsigned int /*zero-based start*/, unsigned char /*region id*/> > &m_regions, 
	const unsigned int &m_loc)
{
	if( m_regions.empty() ){
		return 0;
	}

	vector< pair<unsigned int, unsigned char> >::const_iterator last = m_regions.begin();

	for(vector< pair<unsigned int, unsigned char> >::const_iterator i = m_regions.begin() + 1;i != m_regions.end();++i){

		if(i->first > m_loc){
			return last->second;
		}

		last = i;
	}
	
	return last->second;
}

// Return the kmer Word from m_seq centered at m_center_loc. Return pair.second == true if we extracted a valid
// word and false if we were unable to extract valid word.
const pair<Word, bool> extract_centered_kmer(const unsigned int &m_center_loc, const unsigned int &m_motif_len, 
	const Sequence &m_seq)
{
	pair<Word, bool> ret(0x0, false);

	// Make sure that m_motif_len is odd (since we are extracting the center base)
	if(m_motif_len%2 != 1){
		throw __FILE__ ":extract_centered_kmer: The motif length must be odd";
	}

	// The number of bases on either side of the center base
	const unsigned int num_flank = (m_motif_len - 1)/2;

	if( (m_center_loc < num_flank) || ( (m_center_loc + num_flank) > m_seq.size() ) ){

		// The desired kmer cannot be extracted since it extends past the edge of the sequence
		return ret;
	}

	const unsigned int begin = m_center_loc - num_flank;
	const unsigned int end = m_center_loc + num_flank;

	for(unsigned int i = begin;i <= end;++i){
		
		switch( m_seq[i].base() ){
			case SequenceBase::A:
				ret.first = (ret.first << 2) | BASE_A;
				break;
			case SequenceBase::T:
				ret.first = (ret.first << 2) | BASE_T;
				break;
			case SequenceBase::G:
				ret.first = (ret.first << 2) | BASE_G;
				break;
			case SequenceBase::C:
				ret.first = (ret.first << 2) | BASE_C;
				break;
			default:
				// Invalid base
				return ret;
		};
	}

	// This is a valid word
	ret.second = true;

	return ret;
}

Tensor<float> load_param(const string &m_filename, const TensorDim &m_dim, const size_t &m_num_word)
{
	Tensor<float> ret(m_dim.num_codon_loc, 
				m_dim.num_region, 
				m_dim.num_mutation,
				m_num_word);

	ifstream fin( m_filename.c_str() );

	if(!fin){
		throw __FILE__ ":load_param: Unable to open parameter file";
	}

	string line;

	unsigned int mutation = NO_MUTATION;
	unsigned int codon_loc = 0;
	unsigned int region_loc = 0;

	while( getline(fin, line) ){

		// Is this a header
		if(line.find("k-mer motif for mutation") != string::npos){
			
			// Replace all of the ';' with spaces to make token extraction easier
			for(string::iterator i = line.begin();i != line.end();++i){
				if(*i == ';'){
					*i = ' ';
				}
			}

			const vector<string> header = split(line);

			if(header.size() != 10){
				throw __FILE__ ":load_param: Did not read the expected number of head columns";
			}

			// Column 4 is the mutation type
			size_t mutation_col = 4;
			size_t codon_col = 7;
			size_t region_col = 9;

			if(header[mutation_col] == "->A"){
				mutation = MUTATION_SUB_A;
			}
			else{
				if(header[mutation_col] == "->T"){
					mutation = MUTATION_SUB_T;
				}
				else{
					if(header[mutation_col] == "->G"){
						mutation = MUTATION_SUB_G;
					}
					else{
						if(header[mutation_col] == "->C"){
							mutation = MUTATION_SUB_C;
						}
						else{
							if(header[mutation_col] == "del"){
								mutation = MUTATION_DELETION;
							}
							else{
								if(header[mutation_col] == "ins"){
									mutation = MUTATION_INSERTION;
								}
								else{
									throw __FILE__ ":load_param: Unknown mutation type";
								}
							}
						}
					}
				}
			}

			if(mutation > m_dim.num_mutation){
				throw __FILE__ ":load_param: Invalid mutation";
			}

			codon_loc = abs( atoi( header[codon_col].c_str() ) );

			if(codon_loc >= m_dim.num_codon_loc){
				throw __FILE__ ":load_param: Invalid codon value";
			}

			region_loc = atoi( header[region_col].c_str() );

			if(region_loc >= m_dim.num_region){
				throw __FILE__ ":load_param: Invalid region value";
			}
		}
		else{

			const vector<string> data = split(line);

			if(data.size() != 2){
				continue;
			}

			const size_t word_col = 0;
			const size_t val_col = 1;

			Word w = 0;
			
			for(string::const_iterator i = data[word_col].begin();i != data[word_col].end();++i){

				switch(*i){
					case 'A':
						w = (w << 2) | BASE_A;
						break;
					case 'T':
						w = (w << 2) | BASE_T;
						break;
					case 'G':
						w = (w << 2) | BASE_G;
						break;
					case 'C':
						w = (w << 2) | BASE_C;
						break;
					default:
						throw __FILE__ ":load_param: Invalid base in word";
						break;
				};
			}

			if(w >= m_num_word){
				throw __FILE__ ":load_param: Word out of bounds";
			}

			// Add a random perturbation to test the effect of ties.
			//const float r = 0.1*float(rand())/RAND_MAX; <-- Does not have much effect!!
			//ret(codon_loc, region_loc, mutation, w) = atof(data[ val_col].c_str() ) + r;

			ret(codon_loc, region_loc, mutation, w) = atof(data[ val_col].c_str() );
		}
	}

	return ret;
}