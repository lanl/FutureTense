#ifndef __INDEX
#define __INDEX

#define		MUTATION_SUB_A		0	// == BASE_A in word.h
#define		MUTATION_SUB_C		1	// == BASE_C in word.h
#define		MUTATION_SUB_G		2	// == BASE_G in word.h
#define		MUTATION_SUB_T		3	// == BASE_T in word.h
#define		MUTATION_DELETION	4
#define		MUTATION_INSERTION 	5
#define		NUM_MUTATION		6
#define		NO_MUTATION 		NUM_MUTATION

// The mutation names *must* match the order of the above definition
const char* mutation_name[] = {"->A", "->C", "->G", "->T", "del", "ins", "none"};

// Combine:
// 		1) An index into an array of existing word score. The index is the numeric kmer value.
//		2) A single mutation state
//		3) The region that this word is found in
//		4) The weight == number of observations with matching word, state and region
struct WordIndex
{
	unsigned int weight;

	// The use of a short unsigned int limits the maximum motif length to 8 (i.e. 4^8 = 2^16)
	//short unsigned int index; // <-- the index is the numeric value of the kmer

	// The use of an unsigned int limits the maximum motif length to 16 (i.e. 4^16 = 2^32)
	unsigned int index; // <-- the index is the numeric value of the kmer

	unsigned char region; // The region of the genome (i.e. gene. orf, UTR, ...)

	// "Hey", you ask, why are bit-fields being used for codon_loc and state? Without bit-fields, sizeof(WordIndex) == 12.
	// Using bit-fields, sizeof(WordIndex) == 8 (which should help preserve precious cache). This speeds up execution:
	// 		Bit-field time for one search iteration: approx 40 sec
	//		Full bit time for one search iteration: approx 45 sec
    
	//unsigned char codon_loc : 3; // The codon location of the center of the word, either 1, 2, 3 or 0 (for no codon)
	//unsigned char state : 3; // The mutation (including NO_MUTATION) that is associated with this word & region combination
	//unsigned char center_base : 2;

	// Current benchmarks (Apr 28, 2022) suggest no speed up using bit-fields, so don't add the extract complexity!
	unsigned char codon_loc; // The codon location of the center of the word, either 1, 2, 3 or 0 (for no codon)
	unsigned char state; // The mutation (including NO_MUTATION) that is associated with this word & region combination
	unsigned char center_base;

	WordIndex()
	{
		// Do nothing
	};

	WordIndex(const unsigned int &m_weight, const unsigned short int &m_index,
		const unsigned char &m_region, const unsigned char &m_codon_loc, const unsigned char &m_state,
		const unsigned char &m_center_base) :
		weight(m_weight), index(m_index), region(m_region), codon_loc(m_codon_loc), state(m_state), center_base(m_center_base)
	{
		if(m_weight == 0){
			throw __FILE__ ":WordIndex(): Invalid weight";
		}

		if(m_state > NO_MUTATION){
			throw __FILE__ ":WordIndex(): Invalid mutation state";
		}

		if(m_codon_loc > 3){
			throw __FILE__ ":WordIndex(): Invalid codon location";
		}
	};
};


struct sort_by_mutation
{
	// Order WordIndex structures so that the *smaller* mutation state values are at
	// the begining.
	inline bool operator()(const WordIndex &m_a, const WordIndex &m_b) const
	{
		if(m_a.state == m_b.state){

			if(m_a.region == m_b.region){

				if(m_a.codon_loc == m_b.codon_loc){
					return m_a.index < m_b.index;
				}

				return m_a.codon_loc < m_b.codon_loc;
			}

			return m_a.region < m_b.region;
		}

		return m_a.state < m_b.state;
	};
};

#endif //__INDEX