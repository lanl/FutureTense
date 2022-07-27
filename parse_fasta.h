#ifndef __PARSE_FASTA
#define __PARSE_FASTA

#include <string>
#include <unordered_map>
#include "sequence.h"

// Load fasta records
void load_fasta(std::unordered_map<std::string, Sequence> &m_seq, const std::string &m_filename);

#endif // __PARSE_FASTA