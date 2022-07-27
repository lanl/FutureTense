// Output the difference between two, BLAST-aligned sequences (query and subject) as a CIGAR-like string.
// The difference string encodes the changes need to map the query to the subject. In other words, the
// subject sequence provides the reference frame for the difference string.

#include <stdlib.h>
#include <zlib.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include "quasimodo.h"
#include "word.h"
#include "split.h"

using namespace std;

#define		GENOME_DIFF_VERSION		"0.1, Aug 25, 2021"

#define		MAX_GENOME_SIZE			50000

unsigned char iUPAC_base(char m_base);
bool is_insertion(const string &m_buffer);
bool is_gap(const char m_buffer);
bool consistency_check(const string &m_query_accession, const deque< pair<unsigned int, string> > &m_diff);

int main(int argc, char *argv[])
{
	try{

		// Command line options:
		// [-i <input BLAST mapping file>] (default is stdin)
		// [-o <output file>] (default is stdout)

		const char* options = "i:o:?h";
		//int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;

		string input_filename;
		string output_filename;
		bool print_usage = false;

		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

			switch( opt_code ){
				case 0:
					cerr << "Unknown command line flag!" << endl;
					
					print_usage = true;
					break;
					
				case 'i':
					input_filename = optarg;
					break;
				case 'o':
					output_filename = optarg;
					break;
				case 'h':
				case '?':

					print_usage = true;
					break;
				default:
					cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
					print_usage = true;
					break;
			};
		}

		if(print_usage){

			cerr << "Usage for genome-diff: " << GENOME_DIFF_VERSION << endl;
			cerr << "\t[-i <input BLAST alignment>] (default is stdin)" << endl;
			cerr << "\t[-o <output file>] (default is stdout)" << endl;
			
			return EXIT_SUCCESS;
		}

		ifstream file_in;

		if( !input_filename.empty() ){

			file_in.open( input_filename.c_str() );

			if(!file_in){
				throw "Unable to open input mapping file";
			}
		}

		istream &fin = input_filename.empty() ? cin : file_in;

		ofstream file_out;

		if( !output_filename.empty() ){

			file_out.open( output_filename.c_str() );

			if(!file_out){
				throw "Unable to open output file";
			}
		}

		ostream &fout = output_filename.empty() ? cout : file_out;

		string line;

		string query_accession;
		string subject_accession;
		string query_seq;
		unsigned int query_start;
		unsigned int query_stop;

		// Since there can be repeated elements in query and subject, keep track of the query
		// and subject bases that have already been aligned
		vector<bool> query_mask(MAX_GENOME_SIZE, false);
		vector<bool> subject_mask(MAX_GENOME_SIZE, false);

		deque< pair<unsigned int, string> > diff;

		const string progress = "|/-\\";
		size_t progress_index = 0;

		cerr << "Processing:  ";

		fout << "#query_accession\tsubject_accesion\tdifference" << endl;

		// Read the sequence alignment information from a BLAST-formatted alignment
		while( getline(fin, line) ){

			if( line.find("Query=") != string::npos){

				// Is there an existing difference to write?
				if( !query_accession.empty() ){
					
					// Update the progress
					cerr << '\b' << progress[progress_index];
					progress_index = (progress_index + 1) % progress.size();

					fout << query_accession << '\t' << subject_accession << '\t';

					// Since the differences may be accumulated from disjoint alignments, we need to sort
					// before final output
					sort( diff.begin(), diff.end() );

					consistency_check(query_accession, diff);

					for(deque< pair<unsigned int, string> >::const_iterator i = diff.begin();i != diff.end();++i){
						fout << i->first << i->second;
					}

					fout << endl;
				}

				// Reset the difference string
				diff.clear();

				// Reset the sequence masks
				query_mask.assign(MAX_GENOME_SIZE, false);
				subject_mask.assign(MAX_GENOME_SIZE, false);

				const vector<string> col = split(line);

				if(col.size() < 2){
					throw __FILE__ ":main: Unable to extract query accessopm";
				}

				query_accession = col[1];

				continue;
			}

			if(line.find('>') != string::npos){
				
				// Parse the subject accession
				size_t loc = line.find('>');

				if(loc == string::npos){
					throw __FILE__ ":main: Unable to find '>' for subject accession";
				}

				// Skip the '>'
				++loc;

				const size_t line_len = line.size();

				// Skip any white space between the '>' and the start of the accession
				while( (loc < line_len) && isspace(line[loc]) ){
					++loc;
				}

				subject_accession.clear();

				while( (loc < line_len) && !isspace(line[loc]) ){

					subject_accession.push_back(line[loc]);
					++loc;
				}
			}

			if(line.find("No hits found") != string::npos){

				// There are a few NCBI sequences are a a single base (an 'A') followed by poly-N:
				//		LR991323, LR991330, LR991365, LR963366, LR963146, LR963147, LR963149, OU165929
				cerr << "Warning! No hits found for " << query_accession << endl;

				// Skip this sequence
				query_accession.clear();

				continue;
			}

			if(line.find("Query") == 0){

				stringstream ssin(line);

				if( !(ssin >> line) ){ // Use the "line" to throw away the first column
					throw __FILE__ ":main: Unable to read \"Query\"";
				}

				if( !(ssin >> query_start) ){
					throw __FILE__ ":main: Unable to read query_start";
				}

				if(query_start == 0){ // BLAST coordinates are 1-based
					throw __FILE__ ":main: Invalid query_start";
				}

				// Make the query index 0-based
				--query_start;

				if( !(ssin >> query_seq) ){
					throw __FILE__ ":main: Unable to read query_seq";
				}

				if( !(ssin >> query_stop) ){
					throw __FILE__ ":main: Unable to read query_stop";
				}

				if(query_stop == 0){ // BLAST coordinates are 1-based
					throw __FILE__ ":main: Invalid query_stop";
				}

				// Make the subject (i.e. reference) index 0-based
				--query_stop;

				continue;
			}

			if(line.find("Sbjct") == 0){

				string subject_seq;
				unsigned int subject_start;
				unsigned int subject_stop;

				stringstream ssin(line);

				if( !(ssin >> line) ){ // Use the "line" to throw away the first column
					throw __FILE__ ":main: Unable to read \"Sbjct\"";
				}

				if( !(ssin >> subject_start) ){
					throw __FILE__ ":main: Unable to read subject_start";
				}

				if(subject_start == 0){ // BLAST coordinates are 1-based
					throw __FILE__ ":main: Invalid subject_start";
				}

				// Make the subject (i.e. reference) index 0-based
				--subject_start;

				if( !(ssin >> subject_seq) ){
					throw __FILE__ ":main: Unable to read subject_seq";
				}

				if( subject_seq.size() != query_seq.size() ){
					throw __FILE__ ":main: |subject| != |query|";
				}

				if( !(ssin >> subject_stop) ){
					throw __FILE__ ":main: Unable to read subject_stop";
				}

				if(subject_stop == 0){ // BLAST coordinates are 1-based
					throw __FILE__ ":main: Invalid subject_stop";
				}

				// Make the subject (i.e. reference) index 0-based
				--subject_stop;

				// Only accept alignments that have *both* the query and the subject aligned on the *plus* strand
				if( (query_start > query_stop) || (subject_start > subject_stop) ){
					continue;
				}

				string::const_iterator q = query_seq.begin();
				string::const_iterator s = subject_seq.begin();

				while( q != query_seq.end() ){

					// Make the query base upper case
					const char query_base = toupper(*q);

					if( is_gap(*s) ){

						// Gap in subject == insertion in query

						// Is this the continuation of an existing insertion?
						if( !diff.empty() && 
							(subject_start == diff.back().first) ){
							
							if( !islower(diff.back().second[0]) ){
								throw __FILE__ ":main: Did not find expected insertion";
							}

							diff.back().second += tolower(*q);
						}
						else{

							if( (query_start >= MAX_GENOME_SIZE) || (subject_start >= MAX_GENOME_SIZE) ){
								throw __FILE__ ":main: MAX_GENOME_SIZE exceeded!";
							}

							// Don't incude "secondary" alignment (due to repeated sequences)
							if(!query_mask[query_start] && !subject_mask[subject_start]){

								// This is a new insertion. Use lower-case symbols to indicate
								// sequences inserted into the query (relative to the subject)
								diff.push_back( make_pair( subject_start, string() ) );
								diff.back().second.push_back( tolower(query_base) );
							}
						}
					}
					else{
						
						// Resolve any iUPAC ambiguities when looking for differences between
						// the query and subject bases
						if( (iUPAC_base(query_base) & iUPAC_base(*s) ) == 0){

							if( (query_start >= MAX_GENOME_SIZE) || (subject_start >= MAX_GENOME_SIZE) ){
								throw __FILE__ ":main: MAX_GENOME_SIZE exceeded!";
							}

							// Don't incude "secondary" alignments that may bedue to repeated sequences
							if( !subject_mask[subject_start] ){

								// There can be multiple query gaps at the *same* query location (but with
								// different subject locations).
								if(!query_mask[query_start] || is_gap(*q) ){

									const string query_payload = string(1, query_base);

									diff.push_back( make_pair(subject_start, query_payload) );
								}
							}
						}
					}

					if( (query_start >= MAX_GENOME_SIZE) || (subject_start >= MAX_GENOME_SIZE) ){
						throw __FILE__ ":main: MAX_GENOME_SIZE exceeded!";
					}

					query_mask[query_start] = true;
					subject_mask[subject_start] = true;

					if( !is_gap(*s) ){

						// Only increment the subject location if the query does *not*
						// contain an insertion relative to the subject (which will appear
						// as a gap in the subject).
						++subject_start;
					}

					if( !is_gap(*q) ){

						// Only increment the query location if the subject does *not*
						// contain an insertion relative to the query (which will appear
						// as a gap in the query).
						++query_start;
					}

					++q;
					++s;
				}

				continue;
			}
		}

		// Write the last difference
		if( !query_accession.empty() ){
					
			fout << query_accession << '\t' << subject_accession << '\t';

			// Since the differences may be accumulated from disjoint alignments, we need to sort
			// before final output
			sort( diff.begin(), diff.end() );

			consistency_check(query_accession, diff);

			for(deque< pair<unsigned int, string> >::const_iterator i = diff.begin();i != diff.end();++i){
				fout << i->first << i->second;
			}

			fout << endl;
		}

		cerr << '\b' << "Complete!" << endl;
	}
	catch(const char *error){

		cerr << "Caught the error: " << error << endl;
		return EXIT_FAILURE;
	}
	catch(const exception &error){

		cerr << "Caught the error: " << error.what() << endl;
		return EXIT_FAILURE;
	}
	catch(...){

		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

unsigned char iUPAC_base(char m_base)
{
	enum {
		BASE_GAP = 0,
		BASE_A = 1 << 0,
		BASE_T = 1 << 1,
		BASE_G = 1 << 2,
		BASE_C = 1 << 3,
		BASE_M = (BASE_A | BASE_C),
		BASE_R = (BASE_G | BASE_A),
		BASE_S = (BASE_G | BASE_C),
		BASE_V = (BASE_G | BASE_C | BASE_A),
		BASE_W = (BASE_A | BASE_T),
		BASE_Y = (BASE_T | BASE_C),
		BASE_H = (BASE_A | BASE_C | BASE_T),
		BASE_K = (BASE_G | BASE_T),
		BASE_D = (BASE_G | BASE_A | BASE_T),
		BASE_B = (BASE_G | BASE_T | BASE_C),
		BASE_N = (BASE_A | BASE_T | BASE_C | BASE_G)
	};

	switch(m_base){
		case 'A': case 'a':
			return BASE_A;
		case 'T': case 't':
		case 'U': case 'u':
			return BASE_T;
		case 'G': case 'g':
			return BASE_G;
		case 'C': case 'c':
			return BASE_C;
		case 'M': case 'm':
			return BASE_M;
		case 'R': case 'r':
			return BASE_R;
		case 'S': case 's':
			return BASE_S;
		case 'V': case 'v':
			return BASE_V;
		case 'W': case 'w':
			return BASE_W;
		case 'Y': case 'y':
			return BASE_Y;
		case 'H': case 'h':
			return BASE_H;
		case 'K': case 'k':
			return BASE_K;
		case 'D': case 'd':
			return BASE_D;
		case 'B': case 'b':
			return BASE_B;
		case 'N': case 'n':
			return BASE_N;
		case '-':
			return BASE_GAP;
		default:
			cerr << "Invalid base = " << m_base << " (" << int(m_base) << ")" << endl;
			throw __FILE__ ":iUPAC_base: Illegal base";
			break;
	};

	// We should never get here!
	throw __FILE__ ":iUPAC_base: Invalid base";
	return BASE_GAP;
}

#ifdef USE_PNG
// For libPNG
void png_user_error(png_struct* png_ptr, const char* m_error)
{
	fprintf(stderr, "Error writing PNG image: %s\n", m_error);
};

void png_user_warning(png_struct* png_ptr, const char* m_warning)
{
	fprintf(stderr, "Warning writing PNG image: %s\n", m_warning);
};
#endif // USE_PNG

bool is_insertion(const string &m_buffer)
{
	for(string::const_iterator i = m_buffer.begin();i != m_buffer.end();++i){
		if( (*i < 'a') || (*i > 'z') ){
			return false;
		}
	}

	return true;
}

bool is_gap(const char m_buffer)
{
	return (m_buffer == '-');
}

bool consistency_check(const string &m_query_accession, const deque< pair<unsigned int, string> > &m_diff)
{
	bool is_valid = true;

	// Check for a consistent set of differences. Since there can be multiple, disjoint alignments
	// for a single query genome, are all of the predicted differences to the reference genome
	// self-consistent?
	deque< pair<unsigned int, string> >::const_iterator last_iter = m_diff.begin();

	for(deque< pair<unsigned int, string> >::const_iterator i = m_diff.begin();i != m_diff.end();++i){

		if( (i->first == last_iter->first) && (i->second != last_iter->second) ){

			// Only test if both difference are *not* insertions or *are* insertions
			if( is_insertion(i->second) == is_insertion(last_iter->second) ){

				// Insertions can have the same coordinates as substitutions. For example:
				//	TGTTCT								TGTTCT
				//	||   |	==> is ambiguous with ==>	||   |
				//	TG--AT								TGA--T
				cerr << "Inconsistent alignment for " << m_query_accession << endl;
				cerr << "\tlast: " << last_iter->first << '\t' << last_iter->second << endl;
				cerr << "\tcurr: " << i->first << '\t' << i->second << endl;

				is_valid = false;
			}
		}

		last_iter = i;
	}

	return is_valid;
}	
