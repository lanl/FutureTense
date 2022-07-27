#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "genome_difference.h"

using namespace std;

struct sort_by_difference
{
	inline bool operator()(const pair<string, string> &m_a, const pair<string, string> &m_b) const
	{
		return (m_a.second < m_b.second);	
	};
};

struct unique_by_difference
{
	inline bool operator()(const pair<string, string> &m_a, const pair<string, string> &m_b) const
	{
		return (m_a.second == m_b.second);	
	};
};

struct sort_by_distance_to_reference
{
	inline bool operator()(const pair<string, GenomeDifference> &m_a, const pair<string, GenomeDifference> &m_b) const
	{
		return ( m_a.second.edit_distance_to_reference() < m_b.second.edit_distance_to_reference() );	
	};
};

int main(int argc, char *argv[])
{
	try{

		if(argc != 2){

			cerr << "Usage: " << argv[0] << " <genome diff file>" << endl;
			return EXIT_SUCCESS;
		}

		const string filename = argv[1];

		ifstream fin( filename.c_str() );

		if(!fin){
			cerr << "Unable to open genome diff file, " << filename << ", for reading" << endl;
		}

		// Check to make sure that all genomes are being compared to the same reference genome
		string reference_accession;

		deque< pair<string /*accession*/, string /*difference*/> > data;

		string line;

		while(getline(fin, line)){

			size_t comment_loc = line.find('#');

			if(comment_loc != string::npos){
				line = line.substr(0, comment_loc);
			}

			if( line.empty() ){
				continue;
			}

			stringstream ssin(line);

			string query_accession;
			string subject_accession;
			string buffer;

			if( !(ssin >> query_accession) ){

				cerr << "Unable to read query accession" << endl;
				return EXIT_FAILURE;
			}

			if( !(ssin >> subject_accession) ){

				cerr << "Unable to read subject accession" << endl;
				return EXIT_FAILURE;
			}

			// Make sure that all genome have been compared to the *same* reference genome
			if( reference_accession.empty() ){
				reference_accession = subject_accession;
			}

			if(reference_accession != subject_accession){
				throw __FILE__ ":main: Invalid subject accession; All genomes must use the same reference";
			}

			ssin >> buffer;

			// Entries with empty buffer (i.e. no difference string) as identical to the reference genome
			if( !buffer.empty() ){
				data.push_back( make_pair(query_accession, buffer) );
			}
		}

		cerr << "Found " << data.size() << " genomes that differ from the reference" << endl;

		// Sort the data lexographically by different string
		sort( data.begin(), data.end(), sort_by_difference() );

		data.erase( unique( data.begin(), data.end(), unique_by_difference() ), data.end() );

		cerr << "Found " << data.size() << " unique genomes" << endl;

		// Count the number of unique substitutions, insertions and deletions
		unordered_map<string, unsigned int> substitution_count;
		unordered_map<string, unsigned int> deletetion_count;
		unordered_map<string, unsigned int> insertion_count;

		for(deque< pair<string, string> >::const_iterator i = data.begin();i != data.end();++i){

			string buffer;
			char last = '?';

			for(string::const_iterator j = i->second.begin();j != i->second.end();++j){

				if( ( (last < '0') || (last > '9') ) && (*j >= '0') && (*j <= '9') ){

					if( !buffer.empty() ){
						if(buffer.back() == '-'){
							++deletetion_count[buffer];
						}
						else{
							if( isupper( buffer.back() ) ){
								++substitution_count[buffer];
							}
							else{
								++insertion_count[buffer];
							}
						}
					}

					buffer.clear();
				}

				buffer.push_back(*j);
				last = *j;
			}

			if( !buffer.empty() ){
				if(buffer.back() == '-'){
					++deletetion_count[buffer];
				}
				else{
					if( isupper( buffer.back() ) ){
						++substitution_count[buffer];
					}
					else{
						++insertion_count[buffer];
					}
				}
			}
		}

		cerr << "Found " << substitution_count.size() << " unique substitutions" << endl;
		cerr << "Found " << deletetion_count.size() << " unique deletions" << endl;
		cerr << "Found " << insertion_count.size() << " unique insertions" << endl;

		// Print the substitution counts
		//for(unordered_map<string, unsigned int>::const_iterator i = substitution_count.begin();
		//	i != substitution_count.end();++i){
		//
		//	cout << i->second << '\t' << i->first.substr(0, i->first.size() - 1) << '\t' << i->first.back() << endl;
		//}

		// Find the closest genome that has a smaller edit distance to the reference
		deque< pair<string, GenomeDifference> > diff;

		for(deque< pair<string, string> >::const_iterator i = data.begin();i != data.end();++i){

			diff.push_back( make_pair( i->first, GenomeDifference(i->second) ) );
		}

		// Sort the genomes by their edit distance to the reference
		sort( diff.begin(), diff.end(), sort_by_distance_to_reference() );

		const size_t num_genome = diff.size();

		// Precompute the edit distance from each genome to the reference
		vector<unsigned int> edit_distance(num_genome);

		for(size_t i = 0;i < num_genome;++i){
			edit_distance[i] = diff[i].second.edit_distance_to_reference();
		}

		// Count the number of children beloning to each genome
		unordered_map<string, unsigned int> num_children;

		cerr << "Progress: ";

		const size_t update_every = max(1UL, num_genome/100);
		string info;

		for(size_t i = 0;i < num_genome;++i){
			
			if(i%update_every == 0){

				stringstream ssin;

				ssin << (100.0*i)/num_genome << '%';

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

			const unsigned int edit_distance_i = edit_distance[i];

			size_t parent_index = num_genome;
			unsigned int parent_dist = edit_distance_i;

			for(size_t j = 0;j < num_genome;++j){

				if(edit_distance[j] >= edit_distance_i){
					break;
				}

				const unsigned int d = diff[i].second.distance(diff[j].second);

				if(d < parent_dist){

					parent_dist = d;
					parent_index = j;
				}
			}

			if(parent_index == num_genome){

				cout << "reference\tparent_of\t" << diff[i].first << endl;

				++num_children["reference"];
			}
			else{
				cout << diff[parent_index].first << "\tparent_of\t" << diff[i].first << endl;

				++num_children[diff[parent_index].first];
			}
		}

		cerr << " complete! Summarizing tree" << endl;

		// Sort each genome by its number of children
		deque< pair<unsigned int, string> > child_ranked_genomes;

		for(unordered_map<string, unsigned int>::const_iterator i = num_children.begin();i != num_children.end();++i){
			child_ranked_genomes.push_back( make_pair(i->second, i->first) );
		}

		sort( child_ranked_genomes.begin(), child_ranked_genomes.end() );

		for(deque< pair<unsigned int, string> >::const_reverse_iterator i = child_ranked_genomes.rbegin();
			i != child_ranked_genomes.rend();++i){

			cout << "#\t" << i->second << '\t' << i->first << endl;
		}

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
