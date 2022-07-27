#include "parse_codon.h"
#include <zlib.h>
#include <string.h>
#include <iostream>

using namespace std;

void insert_codon_loc(const string &m_filename, MAP<string /*accession*/, Sequence> &m_db)
{
    // Use zlib to read both compressed and uncompressed codon mapping files.
	gzFile fin = gzopen(m_filename.c_str(), "r");
	
	if(fin == NULL){
		
		cerr << "Error opening: " << m_filename << endl;
		throw __FILE__ ":insert_codon_loc: Unable to open file";
	}

	const unsigned int buffer_len = 65536;
	char *buffer = new char[buffer_len];
	
    // Count the number of sequences for which we were able to successfully apply the codon mapping
    size_t count = 0;

	while( gzgets(fin, buffer, buffer_len) ){
		
        size_t len = strlen(buffer);

        if(len == 0){
            continue;
        }

        if(buffer[len - 1] != '\n'){
            throw __FILE__ ":insert_codon_loc: Line buffer exhausted -> line is too long!";
        }

        // Remove any end of line symbol(s)
        while(len > 0){
            
            const char c = buffer[len - 1];

            if( (c == '\n') || (c == '\r') ){

                buffer[len - 1] = '\0';
                --len;
                continue;
            }

            break;
        }

        // Skip empty lines
        if(len == 0){
            continue;
        }

        string accession;

        // Extract the accession
        char* ptr = buffer;

        // Skip any leading white space
        while( isspace(*ptr) && (*ptr != '\0') ){
            ++ptr;
        }

        // Skip lines that only contain white space
        if(*ptr == '\0'){
            continue;
        }

        while( !isspace(*ptr) && (*ptr != '\0') ){

            accession.push_back(*ptr);
            ++ptr;
        }

        if( accession.empty() ){
            throw __FILE__ ":insert_codon_loc: Unable to extract accession";
        }

        // Look up the accession
        MAP<string /*accession*/, Sequence>::iterator iter = m_db.find(accession);

        if( iter == m_db.end() ){
            
            cerr << "Missing accession: " << accession << endl;
            throw __FILE__ ":insert_codon_loc: Unable to find accession!";
        }

        // Skip the white space that separates the accession from the condon map
        while( isspace(*ptr) && (*ptr != '\0') ){
            ++ptr;
        }

        Sequence::iterator seq_iter = iter->second.begin();

        // Set the per nucleotide codon positions
        while(*ptr != '\0'){

            if( seq_iter == iter->second.end() ){
                throw __FILE__ ":insert_codon_loc: Sequence buffer overflow!";
            }

            switch(*ptr){
                case '0':
                case '1':
                case '2':
                case '3':
                    seq_iter->set_codon(*ptr - '0');
                    break;
                default:
                    throw __FILE__ ":insert_codon_loc: Invalid codon location symbol";
            };

            ++ptr;
            ++seq_iter;
        }

        // Make sure we set the codon location for all bases
        if( seq_iter != iter->second.end() ){
            throw __FILE__ ":insert_codon_loc: Codon mapping underflow";
        }

        ++count;

    }
	
    delete [] buffer;
	gzclose(fin);

    if( count != m_db.size() ){
        throw __FILE__ ":insert_codon_loc: Failed to map codon locations for all database sequences";
    }
}