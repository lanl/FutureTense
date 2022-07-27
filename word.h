#ifndef __KMER_TOOLBOX
#define __KMER_TOOLBOX

#include <string>

// Storing words in 64-bits (at two bits per base) yeilds a max word
// size of 32.
#define		MAX_WORD_LEN	32
typedef		size_t	Word;

// A two bit binary base encoding. Order the bases lexographically
enum {BASE_A, BASE_C, BASE_G, BASE_T, NUM_BASE};

// The base_name *must* match the order of the enumerated bases above
const char base_name[] = "ACGT"; // Bases in lexographic order

////////////////////////////////////////////////////////////////////////////////////
// In word.cpp
std::string word_to_string(const Word &m_w, const size_t &m_k);
Word kmer_word_mask(const size_t &m_len);

// Macros for digesting sequences into kmers
#define ForEachDuplexWord(__SEQ, __LEN)\
	{\
		const Word __comp_shift = 2*(__LEN - 1);\
		const Word __mask = kmer_word_mask(__LEN);\
		Word __w = 0;\
		Word __comp_w = 0;\
		unsigned int __word_len = 0;\
		size_t __k = __LEN;\
		size_t __index = 0;\
		for(Sequence::const_iterator __i = __SEQ.begin();__i != __SEQ.end();++__i,++__index){ \
			++__word_len;\
			switch(__i->base()){\
				case SequenceBase::A:\
					__w = (__w << 2) | BASE_A;\
					__comp_w = (__comp_w >> 2) | (Word(BASE_T) << __comp_shift);\
					break;\
				case SequenceBase::T:\
					__w = (__w << 2) | BASE_T;\
                    __comp_w = (__comp_w >> 2) | (Word(BASE_A) << __comp_shift);\
					break;\
				case SequenceBase::G:\
					__w = (__w << 2) | BASE_G;\
                    __comp_w = (__comp_w >> 2) | (Word(BASE_C) << __comp_shift);\
					break;\
				case SequenceBase::C:\
					__w = (__w << 2) | BASE_C;\
                    __comp_w = (__comp_w >> 2) | (Word(BASE_G) << __comp_shift);\
					break;\
				default:\
					__word_len = 0;\
					break;\
			};

#define ForEachSenseWord(__SEQ, __LEN)\
	{\
		const Word __mask = kmer_word_mask(__LEN);\
		Word __w = 0;\
		unsigned int __word_len = 0;\
		size_t __k = __LEN;\
		size_t __index = 0;\
		for(Sequence::const_iterator __i = __SEQ.begin();__i != __SEQ.end();++__i,++__index){ \
			++__word_len;\
			switch(__i->base()){\
				case SequenceBase::A:\
					__w = (__w << 2) | BASE_A;\
					break;\
				case SequenceBase::T:\
					__w = (__w << 2) | BASE_T;\
					break;\
				case SequenceBase::G:\
					__w = (__w << 2) | BASE_G;\
					break;\
				case SequenceBase::C:\
					__w = (__w << 2) | BASE_C;\
					break;\
				default:\
					__word_len = 0;\
					break;\
			};

#define ForEachAntisenseWord(__SEQ, __LEN)\
	{\
		const Word __comp_shift = 2*(__LEN - 1);\
		const Word __mask = kmer_word_mask(__LEN);\
		Word __comp_w = 0;\
		unsigned int __word_len = 0;\
		size_t __k = __LEN;\
		size_t __index = 0;\
		for(Sequence::const_iterator __i = __SEQ.begin();__i != __SEQ.end();++__i,++__index){ \
			++__word_len;\
			switch(__i->base()){\
				case SequenceBase::A:\
					__comp_w = (__comp_w >> 2) | (Word(BASE_T) << __comp_shift);\
					break;\
				case SequenceBase::T:\
                    __comp_w = (__comp_w >> 2) | (Word(BASE_A) << __comp_shift);\
					break;\
				case SequenceBase::G:\
                    __comp_w = (__comp_w >> 2) | (Word(BASE_C) << __comp_shift);\
					break;\
				case SequenceBase::C:\
                    __comp_w = (__comp_w >> 2) | (Word(BASE_G) << __comp_shift);\
					break;\
				default:\
					__word_len = 0;\
					break;\
			};
				
#define	ErrorWord		(__word_len == 0)
#define	ValidWord		(__word_len >= __k)
#define	SenseWord		(__w & __mask)
#define	AntisenseWord	(__comp_w & __mask)
#define	CanonicalWord	std::min(SenseWord, AntisenseWord)
#define	CurrentBase		(*__i)
#define	Loc3			(__index)				// Location of the *last* base in current kmer
#define	Loc5			( (__index + 1) - __k)	// Location of the *first* base in the current kmer

#define EndWord\
		}\
	}

#endif // __KMER_TOOLBOX
