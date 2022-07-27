#include "parse_fasta.h"
#include <zlib.h>
#include <string.h>
#include <iostream>
#include <deque>

using namespace std;

string extract_accession(const string &m_defline);

void load_fasta(unordered_map<string, Sequence> &m_db, const std::string &m_filename)
{
	// Use zlib to read both compressed and uncompressed fasta files.
	gzFile fin = gzopen(m_filename.c_str(), "r");
	
	if(fin == NULL){
		
		cerr << "Error opening: " << m_filename << endl;
		throw __FILE__ ":load_fasta: Unable to open fasta file";
	}

	const int buffer_len = 2048;
	char buffer[buffer_len];
	
	string accession;
	Sequence* seq_ptr = NULL;

	while( gzgets(fin, buffer, buffer_len) ){
		
		if(strchr(buffer, '>') != NULL){

			// Remove any end of line symbols
			for(char* p = buffer;*p != '\0';++p){
				if( (*p == '\n') || (*p == '\r') ){
					*p = '\0';
				}
			}

			accession = extract_accession(buffer);

			// Check for duplicates
			if( m_db.find(accession) != m_db.end() ){

				cerr << "Duplicate accession: " << accession << endl;
				throw __FILE__ ":load_fasta: Found a duplicate accession";
			}

			m_db[accession] = Sequence();

			seq_ptr = &(m_db[accession]);
		}
		else{

			if(seq_ptr == NULL){
				continue;
			}

			for(char* p = buffer;*p != '\0';++p){
			
				// Original version:
				// Skip spaces and gap symbols
				//if( !isspace(*p) && !(*p == '-') ){
					//seq_ptr->push_back( SequenceBase(*p) );
				//}

				// Updated version to skip some unneeded tests and initializations
				switch(*p){
					case 'A': case 'a':
						seq_ptr->push_back( SequenceBase(SequenceBase::A) );
						break;
					case 'C': case 'c':
						seq_ptr->push_back( SequenceBase(SequenceBase::C) );
						break;
					case 'G': case 'g':
						seq_ptr->push_back( SequenceBase(SequenceBase::G) );
						break;
					case 'T': case 't':
						seq_ptr->push_back( SequenceBase(SequenceBase::T) );
						break;
					case 'M': case 'm':
						seq_ptr->push_back( SequenceBase(SequenceBase::M) );
						break;
					case 'R': case 'r':
						seq_ptr->push_back( SequenceBase(SequenceBase::R) );
						break;
					case 'S': case 's':
						seq_ptr->push_back( SequenceBase(SequenceBase::S) );
						break;
					case 'V': case 'v':
						seq_ptr->push_back( SequenceBase(SequenceBase::V) );
						break;
					case 'W': case 'w':
						seq_ptr->push_back( SequenceBase(SequenceBase::W) );
						break;
					case 'Y': case 'y':
						seq_ptr->push_back( SequenceBase(SequenceBase::Y) );
						break;
					case 'H': case 'h':
						seq_ptr->push_back( SequenceBase(SequenceBase::H) );
						break;
					case 'K': case 'k':
						seq_ptr->push_back( SequenceBase(SequenceBase::K) );
						break;
					case 'D': case 'd':
						seq_ptr->push_back( SequenceBase(SequenceBase::D) );
						break;
					case 'B': case 'b':
						seq_ptr->push_back( SequenceBase(SequenceBase::B) );
						break;
					case 'N': case 'n':
						seq_ptr->push_back( SequenceBase(SequenceBase::N) );
						break;
				};
			}
		}
	}
	
	gzclose(fin);
}

string extract_accession(const string &m_defline)
{
	size_t len = m_defline.size();

	size_t start = m_defline.find('>');

	if(start == string::npos){
		throw __FILE__ ":extract_accession: Did not find '>'";
	}

	// Skip the '>'
	++start;

	while( (start < len) && isspace(m_defline[start]) ){
		++start;
	}

	size_t stop = start;

	while( (stop < len) && !isspace(m_defline[stop])  ){
		++stop;
	}

	return m_defline.substr(start, stop - start);
}