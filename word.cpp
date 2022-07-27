#include "word.h"

#include <deque>
#include <algorithm>

using namespace std;

Word kmer_word_mask(const size_t &m_len)
{
	Word ret = 0;
	
	// The following code is commented out because it does not
	// work when the m_kmer_len is equal to full word length
	// (we overflow the Word).
        //const Word ret = ( 1UL << (2*m_len) ) - 1;

	for(size_t i = 0;i < 2*m_len;++i){
		ret |= (1UL << i);
	}

	return ret;
}

string word_to_string(const Word &m_w, const size_t &m_k)
{
	Word w = m_w;
	deque<char> buffer;
	
	for(size_t i = 0;i < m_k;++i){
		
		switch(w & 3){
			case BASE_A:
				buffer.push_front('A');
				break;
			case BASE_T:
				buffer.push_front('T');
				break;
			case BASE_G:
				buffer.push_front('G');
				break;
			case BASE_C:
				buffer.push_front('C');
				break;
			default:
				throw __FILE__ ":word_to_string: Unknown base!";
		};
		
		w = (w >> 2);
	}
	
	return string( buffer.begin(), buffer.end() );
}
