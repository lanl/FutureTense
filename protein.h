#ifndef __PROTEIN_SEQ
#define __PROTEIN_SEQ

#include <vector>
#include <string>

std::vector< std::pair<char /*amino acid*/, unsigned char /*codon position*/> > protein_seq(const std::string &m_rna,
	const size_t &m_min_len);

#endif // __PROTEIN_SEQ