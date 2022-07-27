#ifndef __PARSE_CODON
#define __PARSE_CODON

#include <string>
#include "map.h"
#include "sequence.h"

void insert_codon_loc(const std::string &m_filename, MAP<std::string /*accession*/, Sequence> &m_db);

#endif // __PARSE_CODON