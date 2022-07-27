#include "protein.h"
#include <iostream> // DEBUG
#include <deque>

using namespace std;

char rna_to_aa(const unsigned int &m_codon);

#define		STOP_CODON		'-'

vector< pair<char /*amino acid*/, unsigned char /*codon position*/> > protein_seq(const string &m_rna, 
	const size_t &m_min_len)
{
	const size_t len = m_rna.size();
	const size_t codon_len = 3;

	vector< pair<char /*amino acid*/, unsigned char /*codon position*/> > ret( len, make_pair(STOP_CODON, 0) );

	// Iteratively find the longest continguous protein sequences
	while(true){

		deque< pair<size_t /*RNA index*/, char /*amino acid*/> > best;
		deque< pair<size_t /*RNA index*/, char /*amino acid*/> > curr;

		for(size_t frame = 0;frame < codon_len;++frame){

			for(size_t i = frame;i < (len - codon_len);i += codon_len){

				if(ret[i].first != STOP_CODON){

					// Previously assigned protein sequence
					if( curr.size() > best.size() ){
						best = curr;
					}

					curr.clear();

					continue;
				}

				unsigned int codon = toupper(m_rna[i]);
				codon = (codon << 8) | toupper(m_rna[i + 1]);
				codon = (codon << 8) | toupper(m_rna[i + 2]);

				const char aa = rna_to_aa(codon);

				if(aa == STOP_CODON){

					// Previously assigned protein sequence
					if( curr.size() > best.size() ){
						best = curr;
					}

					curr.clear();

					continue;
				}

				// Require proteins to start with methoinine
				if(curr.empty() && (aa != 'M') ){
					continue;
				}

				curr.push_back( make_pair(i, aa) );
			}

			curr.clear();
		}

		if(best.size() < m_min_len){
			break;
		}

		// DEBUG
		//cerr << "|Best seq| = " << best.size() << " from " << best.front().first << " to " << best.back().first << endl;

		for(deque< pair<size_t /*RNA index*/, char /*amino acid*/> >::const_iterator i = best.begin();i != best.end();++i){

			if(ret[i->first].first != STOP_CODON){
				throw __FILE__ ":protein_seq: Overwriting previous protein!";
			}

			ret[i->first] = make_pair(i->second, 1);
			ret[i->first + 1] = make_pair(i->second, 2);
			ret[i->first + 2] = make_pair(i->second, 3);
		}
	}

	// DEBUG
	//for(size_t i = 0;i < len;++i){
	//	cout << i << '\t' << m_rna[i] << '\t' << ret[i].first << '\t' << int(ret[i].second) << endl;
	//}

	return ret;
}

char rna_to_aa(const unsigned int &m_codon)
{

	const unsigned int codon_mask = 0x00FFFFFF;

	switch(m_codon & codon_mask){

		case 'TTT': case 'TTC':
			return 'F';

		case 'TTA': case 'TTG': case 'CTT': case 'CTC': case 'CTA': case 'CTG':
			return 'L';
	
		case 'ATT': case 'ATC': case 'ATA':
			return 'I';

		case 'ATG':
			return 'M';

		case 'GTT': case 'GTC': case 'GTA': case 'GTG':
			return 'V';

		case 'TCT': case 'TCC': case 'TCA': case 'TCG':
			return 'S';
		
		case 'CCT': case 'CCC': case 'CCA': case 'CCG':
			return 'P';
	
		case 'ACT': case 'ACC': case 'ACA': case 'ACG':
			return 'T';

		case 'GCT': case 'GCC': case 'GCA': case 'GCG':
			return 'A';
	
		case 'TAT': case 'TAC':
			return 'Y';
	
		case 'TAA': case 'TAG': case 'TGA':
			return STOP_CODON;

		case 'CAT': case 'CAC':
			return 'H';

		case 'CAA': case 'CAG':
			return 'Q';
	
		case 'AAT': case 'AAC':
			return 'N';

		case 'AAA': case 'AAG':
			return 'K';

		case 'GAT': case 'GAC':
			return 'D';
	
		case 'GAA': case 'GAG':
			return 'E';

		case 'TGT': case 'TGC':
			return 'C';
	
		case 'TGG':
			return 'W';
	
		case 'CGT': case 'CGC': case 'CGA': case 'CGG':
			return 'R';
	
		case 'AGT': case 'AGC':
			return 'S';

		case 'AGA': case 'AGG':
			return 'R';

		case 'GGT': case 'GGC': case 'GGA': case 'GGG':
			return 'G';
		
		// There is no default code, since we may encounter degenerate nucleotides that do
		// not have a defined codon (return a STOP_CODON for these).
		//default:
		//
		//	cerr << "Unknown codon: " 
		//		<< char( (m_codon >> 16) & 0xFF) 
		//		<< char( (m_codon >> 8) & 0xFF) 
		//		<< char(m_codon & 0xFF) << endl;
		//
		//	throw __FILE__ ":rna_to_aa: Unknown codon!";
	}

	return STOP_CODON;
}
