#ifndef __SEQ_OVERLAP
#define __SEQ_OVERLAP

#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <vector>
#include <deque>

/////////////////////////////////////////////////////////////////////////////////////////////
//
// SSE instruction defines, macros and datatypes. Note that using SSE4.1 requires
// g++ version >= 4.4 and the compiler flag -msse4.1. Currenlty (3/24/10) the Intel AVX instruction
// set (which offers 256 bit registers) does not apper to be supported in available hardware.
//
// Note that ATA32 requires the sse4.1 instruction set.
// Note that ATA16 requires the sse2 instruction set.
//
// 32-bit OSes require allocated memory be aligned to 16 byte boundaries. This requires use of
// the _mm_malloc/_mm_free functions (instead of new and delete). This dependance does not
// appear to be necessary for 64-bit OSes.
/////////////////////////////////////////////////////////////////////////////////////////////

#include <immintrin.h>

namespace SO {

/////////////////////////////////////////////////////////////////////////////////////////////
// SSE & AVX instruction defines, macros and datatypes. Note that using SSE4.1 requires
// g++ version >= 4.4 and the compiler flag -msse4.1. Using AVX2 requires the compiler 
// flag -mavx2.
/////////////////////////////////////////////////////////////////////////////////////////////
#include <immintrin.h>
	
#ifdef ATA32x8

	// 32-bit OSes require allocated memory be aligned to 16 byte boundaries. This requires use of
	// the _mm_malloc/_mm_free functions (instead of new and delete). This dependance does not
	// appear to be necessary for 64-bit OSes (but is important for backwards compatibility).
	#define	MEMORY_ALIGNMENT			32

	// Use a 32 bit integer data type for dynamic programming. Since AVX provides
	// a 256-bit integer data type, we can process 8 elements in parallel.
	#define		SO_LEN					8u
	#define		SIMD_DATA				__m256i
	
	#define		SIMD_SET1				_mm256_set1_epi32
	#define		SIMD_ADD				_mm256_add_epi32
	#define		SIMD_SUB				_mm256_sub_epi32
	#define		SIMD_CMPEQ				_mm256_cmpeq_epi32
	#define		SIMD_CMPGT				_mm256_cmpgt_epi32
	#define		SIMD_MAX				_mm256_max_epi32
	
	#define		SIMD_BITWISE_AND		_mm256_and_si256
	#define		SIMD_BITWISE_OR			_mm256_or_si256
	#define		SIMD_BITWISE_AND_NOT	_mm256_andnot_si256
	
	#define		SIMD_MALLOC			_mm_malloc
	#define		SIMD_FREE			_mm_free
	
	typedef 	int					SO_Score;
	
#elif ATA32x4
	
	// 32-bit OSes require allocated memory be aligned to 16 byte boundaries. This requires use of
	// the _mm_malloc/_mm_free functions (instead of new and delete). This dependance does not
	// appear to be necessary for 64-bit OSes (but is important for backwards compatibility).
	#define	MEMORY_ALIGNMENT		16

	// Use a 32 bit integer data type for dynamic programming. Since SSE provides
	// a 128-bit integer data type, we can process 4 elements in parallel.
	#define		SO_LEN					4u
	#define		SIMD_DATA				__m128i
	
	#define		SIMD_SET1				_mm_set1_epi32
	#define		SIMD_ADD				_mm_add_epi32
	#define		SIMD_SUB				_mm_sub_epi32
	#define		SIMD_CMPEQ				_mm_cmpeq_epi32
	#define		SIMD_CMPGT				_mm_cmpgt_epi32
	#define		SIMD_MAX				_mm_max_epi32
	
	#define		SIMD_BITWISE_AND		_mm_and_si128
	#define		SIMD_BITWISE_OR			_mm_or_si128
	#define		SIMD_BITWISE_AND_NOT	_mm_andnot_si128
	
	#define		SIMD_MALLOC				_mm_malloc
	#define		SIMD_FREE				_mm_free
	
	typedef 	int					SO_Score;
#elif ATA32x2
	
	// 32-bit OSes require allocated memory be aligned to 16 byte boundaries. This requires use of
	// the _mm_malloc/_mm_free functions (instead of new and delete). This dependance does not
	// appear to be necessary for 64-bit OSes (but is important for backwards compatibility).
	#define	MEMORY_ALIGNMENT			16

	// Use a 32 bit integer data type for dynamic programming. Since MMX provides
	// a 64-bit integer data type, we can process 2 elements in parallel.
	#define		SO_LEN					2u
	#define		SIMD_DATA				__m64
	
	#define		SIMD_SET1				_mm_set1_pi32
	#define		SIMD_ADD				_mm_add_pi32
	#define		SIMD_SUB				_mm_sub_pi32
	#define		SIMD_CMPEQ				_mm_cmpeq_pi32
	#define		SIMD_CMPGT				_mm_cmpgt_pi32
	
	// Note that the MMX instruction set does not provide a built-in max function 
	// We need to provide this functionality ourselves (see _mm_max_pi32 defined below).
	#define		SIMD_MAX				_mm_max_pi32
	
	#define		SIMD_BITWISE_AND		_mm_and_si64
	#define		SIMD_BITWISE_OR			_mm_or_si64
	#define		SIMD_BITWISE_AND_NOT	_mm_andnot_si64
	
	#define		SIMD_MALLOC				_mm_malloc
	#define		SIMD_FREE				_mm_free
	
	typedef 	int						SO_Score;
	
	inline __m64 _mm_max_pi32(__m64 A, __m64 B)
	{
		__m64 mask = _mm_cmpgt_pi32(A, B);
		
		return _mm_or_si64( _mm_and_si64(mask, A), _mm_andnot_si64(mask, B) );
		
	};
	
#elif ATA32x1 // No SIMD!
	
	#define		MEMORY_ALIGNMENT				0

	#define		SO_LEN							1u
	#define		SIMD_DATA						int
	
	#define		SIMD_SET1(_X)					(_X)
	#define		SIMD_ADD(_X, _Y)				( (_X) + (_Y) )
	#define		SIMD_SUB(_X, _Y)				( (_X) - (_Y) )
	#define		SIMD_CMPEQ(_X, _Y)				( ( (_X) == (_Y) ) ? 0xFFFFFFFF : 0x00000000 )
	#define		SIMD_CMPGT(_X, _Y)				( ( (_X) > (_Y) ) ? 0xFFFFFFFF : 0x00000000 )
	
	#define		SIMD_MAX(_X, _Y)				( ((_X) > (_Y)) ? (_X) : (_Y) )
	
	#define		SIMD_BITWISE_AND(_X, _Y)		( (_X) & (_Y) )
	#define		SIMD_BITWISE_OR(_X, _Y)			( (_X) | (_Y) )
	#define		SIMD_BITWISE_AND_NOT(_X, _Y)	( ~(_X) & (_Y) )
	
	#define		SIMD_MALLOC(_X,_Y)				malloc(_X)
	#define		SIMD_FREE						free
	
	typedef 	int								SO_Score;
	
#elif ATA16
	
	// 32-bit OSes require allocated memory be aligned to 16 byte boundaries. This requires use of
	// the _mm_malloc/_mm_free functions (instead of new and delete). This dependance does not
	// appear to be necessary for 64-bit OSes (but is important for backwards compatibility).
	#define	MEMORY_ALIGNMENT					16
	
	// Use a 16 bit integer data type for dynamic programming. Since SSE provides
	// a 128bit integer data type, we can process 8 elements in parallel.
	#define		SO_LEN							8u
	#define		SIMD_DATA						__m128i
	
	#define		SIMD_SET1						_mm_set1_epi16
	#define		SIMD_ADD						_mm_add_epi16
	#define		SIMD_SUB						_mm_sub_epi16
	#define		SIMD_CMPEQ						_mm_cmpeq_epi16
	#define		SIMD_CMPGT						_mm_cmpgt_epi16
	#define		SIMD_MAX						_mm_max_epi16
	
	#define		SIMD_BITWISE_AND				_mm_and_si128
	#define		SIMD_BITWISE_OR					_mm_or_si128
	#define		SIMD_BITWISE_AND_NOT			_mm_andnot_si128
	
	#define		SIMD_MALLOC						_mm_malloc
	#define		SIMD_FREE						_mm_free
	
	typedef 	short int						SO_Score;
#else
	#error "Please specificy a SIMD data size!"
#endif

union sse_elem {
	
	sse_elem()
	{
		// Do nothing
	};
	
	sse_elem(const SO_Score &m_scalar)
	{
		sse = SIMD_SET1(m_scalar);
	};
	
	// If ATA32x1 is defined, then SO_Score == SIMD_DATA. Protect the following
	// constructor to avoid duplication in this case.
	#ifndef ATA32x1
	sse_elem(const SIMD_DATA &m_vector)
	{
		sse = m_vector;
	};
	#endif // ATA32x1

	~sse_elem()
	{
		// Do nothing
	};
	
	SIMD_DATA sse;
	SO_Score v[SO_LEN];
};

//#define		MAX_SEQUENCE_LENGTH	32768
const unsigned int	MAX_SEQUENCE_LENGTH = 32768;

#define		MINIMUM_ALLOWED_SCORE	-9999

// The definition of the amino acids comes *before* the definition of the nucleic acids
// to allow the nucleic acid gap value to be the same as the amino acid gap value.
namespace AA {
	
	// Allowed amino acid bases
	// The actual values must match the score matrix arrangement
	typedef enum {
		A = 0, // Alanine, Ala
		R, // Arginine, Arg
		N, // Asparagine, Asn
		D, // Aspartic acid, Asp
		C, // Cysteine, Cys
		Q, // Glutamine, Gln
		E, // Glutamic acid, Glu
		G, // Glycine, Gly
		H, // Histidine, His
		I, // Isoleucine, Ile
		L, // Leucine, Leu
		K, // Lysine, Lys
		M, // Methionine, Met
		F, // Phenylalanine, Phe
		P, // Proline, Pro
		S, // Serine, Ser
		T, // Threonine, Thr
		W, // Tryptophan, Trp
		Y, // Tyrosine, Tyr
		V, // Valine, Val
		B, // Aspartic acid or Asparagine, Asx
		Z, // Glutamine or Glutamic acid, Glx
		X, // Any amino acid, Xaa
		GAP = (1 << 5) // Make sure that no bits overlap with any nucleic acid

	} amino_acid;
	
	const size_t NUM_AA = (X + 1); // (20 aa + B, Z and X)
	const size_t MATRIX_SIZE = NUM_AA*NUM_AA; // (20 aa + B, Z and X) by (20 aa + B, Z and X)
}

namespace NA {
	
	// Allowed nucleic acid bases
	//	A = adenine
	//	C = cytosine
	//	G = guanine
	//	T = thymine
	//	M = A or C
	//	R = G or A
	//	S = G or C
	//	V = G or C or A
	//	W = A or T
	//	Y = T or C
	//	H = A or C or T
	//	K = G or T
	//	D = G or A or T
	//	B = G or T or C
	//	N = A or T or C or G
	typedef enum {
		A = (1 << 0), 
		C = (1 << 1), 
		G = (1 << 2), 
		T = (1 << 3), 
		M = (A | C),
		R = (G | A),
		S = (G | C),
		V = (G | C | A),
		W = (A | T),
		Y = (T | C),
		H = (A | C | T),
		K = (G | T),
		D = (G | A | T),
		B = (G | T | C),
		N = (A | T | C | G),
		// It make life easier to share the same value for a GAP between
		// nucleic acids and amino acids
		GAP = AA::GAP
	} nucleic_acid;
}

// Both base types (nucleic acid and amino acid) must fit
// into a variable of size base_type
typedef unsigned char base_type;

base_type na_to_bits(const char &m_base);
base_type aa_to_bits(const char &m_base);

class SeqOverlap{

	public:
	
	typedef enum {SmithWaterman, Overlap} AlignmentMode;
	
	private:
		
		struct SO_Elem{

			SO_Elem()
			{
				// Do nothing
			};

			~SO_Elem()
			{
				// Do nothing
			};
			
			// Use the notation of "Biological sequence analysis" by Durbin, Eddy, Krogh and Mitchison
			sse_elem M;
			sse_elem I_query;  // insertion in query
			sse_elem I_target; // insertion in target
									
			// The actual alignment is not needed, just the score and the
			// coordinates of the begining and end of the alignment. For this purpose, propagate
			// the alignment start during the dynamic programming.
			
			// The starting row and column of the alignment
			sse_elem M_start_i;
			sse_elem M_start_j;
			
			sse_elem I_query_start_i;
			sse_elem I_query_start_j;
			
			sse_elem I_target_start_i;
			sse_elem I_target_start_j;
		};
		
		AlignmentMode mode;
		
		SO_Elem *last_row;
		SO_Elem *curr_row;
		SO_Elem max_elem;
		size_t dp_num_elem;
		
		// The location of the max element
		sse_elem stop_i;
		sse_elem stop_j;
		
		sse_elem match;
		sse_elem mask;
		sse_elem mismatch;
		sse_elem gap_existance;
		sse_elem gap_extension;

		// The input sequences in 5'-3' orientation
		sse_elem query[MAX_SEQUENCE_LENGTH];
		sse_elem target[MAX_SEQUENCE_LENGTH];
		
		sse_elem query_len;
		sse_elem target_len;
		
		// The maximum query and target length for the
		// currently packed sequences.
		SO_Score max_query_len;
		SO_Score max_target_len;
		
		bool is_na; // If is_na == false, then we are using amino acids
		
		// If mask_N_na == true, then treat 'N' as a masking character that gets a score of mask
		bool mask_N_na;
		
		std::vector<SO_Score> aa_score_matrix;
		
		void init_BLOSUM62(std::vector<SO_Score> &m_matrix);
		
		inline char bits_to_na(const base_type &m_bits) const
		{
			switch(m_bits){
				case NA::A:
					return 'A';
				case NA::C:
					return 'C';
				case NA::G:
					return 'G';
				case NA::T:
					return 'T';
				case NA::M:
					return 'M';
				case NA::R:
					return 'R';
				case NA::S:
					return 'S';
				case NA::V:
					return 'V';
				case NA::W:
					return 'W';
				case NA::Y:
					return 'Y';
				case NA::H:
					return 'H';
				case NA::K:
					return 'K';
				case NA::D:
					return 'D';
				case NA::B:
					return 'B';
				case NA::N:
					return 'N';
				case NA::GAP:
					return '-';
			};

			throw __FILE__ ":bit_to_na: Unknown base!";
			return 'X'; // Keep the compiler happy
		};
		
		inline char bits_to_aa(const base_type &m_bits) const
		{
			switch(m_bits){
				case AA::A:
					return 'A';
				case AA::R:
					return 'R';
				case AA::N:
					return 'N';
				case AA::D:
					return 'D';
				case AA::C:
					return 'C';
				case AA::Q:
					return 'Q';
				case AA::E:
					return 'E';
				case AA::G:
					return 'G';
				case AA::H:
					return 'H';
				case AA::I:
					return 'I';
				case AA::L:
					return 'L';
				case AA::K:
					return 'K';
				case AA::M:
					return 'M';
				case AA::F:
					return 'F';
				case AA::P:
					return 'P';
				case AA::S:
					return 'S';
				case AA::T:
					return 'T';
				case AA::W:
					return 'W';
				case AA::Y:
					return 'Y';
				case AA::V:
					return 'V';
				case AA::B:
					return 'B';
				case AA::Z:
					return 'Z';
				case AA::X:
					return 'X';
				case AA::GAP:
					return '-';
			};

			throw __FILE__ ":bit_to_aa: Unknown base!";
			return '?'; // Keep the compiler happy
		};

		void translate_na_seq(const std::string &m_input, std::vector<base_type> &m_output)
		{
			#ifdef _DEBUG
			if(!is_na){
				throw __FILE__ ":translate_na_seq: Only valid for NA sequences";
			}
			
			if( m_input.size() != m_output.size() ){
				throw __FILE__ ":translate_na_seq: input/output size mismatch";
			}
			#endif // _DEBUG
			
			std::string::const_iterator i = m_input.begin();
			std::vector<base_type>::iterator o = m_output.begin();

			for(;i != m_input.end();i++, o++){
			
				switch(*i){
					case 'A': case 'a':
						*o = NA::A;
						break;
					case 'T': case 't':
						*o = NA::T;
						break;
					case 'G': case 'g':
						*o = NA::G;
						break;
					case 'C': case 'c':
						*o = NA::C;
						break;
					case 'M': case 'm':
						*o = NA::M;
						break;
					case 'R': case 'r':
						*o = NA::R;
						break;
					case 'S': case 's':
						*o = NA::S;
						break;
					case 'V': case 'v':
						*o = NA::V;
						break;
					case 'W': case 'w':
						*o = NA::W;
						break;
					case 'Y': case 'y':
						*o = NA::Y;
						break;
					case 'H': case 'h':
						*o = NA::H;
						break;
					case 'K': case 'k':
						*o = NA::K;
						break;
					case 'D': case 'd':
						*o = NA::D;
						break;
					case 'B': case 'b':
						*o = NA::B;
						break;
					case 'N': case 'n':
					case 'I': case 'i': // For now, treat inosine as an 'N'
						*o = NA::N;
						break;
					default:
						throw __FILE__ ":translate_na_seq: Illegal base";
						break;
				};
			}
		};

		void translate_aa_seq(const std::string &m_input, std::vector<base_type> &m_output)
		{
			#ifdef _DEBUG
			if(is_na){
				throw __FILE__ ":translate_aa_seq: Only valid for AA sequences";
			}
			
			if( m_input.size() != m_output.size() ){
				throw __FILE__ ":translate_aa_seq: input/output size mismatch";
			}
			#endif // _DEBUG
			
			std::string::const_iterator i = m_input.begin();
			std::vector<base_type>::iterator o = m_output.begin();

			for(;i != m_input.end();i++, o++){
			
				switch(*i){
					case 'A': case 'a':
						*o = AA::A;
						break;
					case 'R': case 'r':
						*o = AA::R;
						break;
					case 'N': case 'n':
						*o = AA::N;
						break;
					case 'D': case 'd':
						*o = AA::D;
						break;
					case 'C': case 'c':
						*o = AA::C;
						break;
					case 'Q': case 'q':
						*o = AA::Q;
						break;
					case 'E': case 'e':
						*o = AA::E;
						break;
					case 'G': case 'g':
						*o = AA::G;
						break;
					case 'H': case 'h':
						*o = AA::H;
						break;
					// J = I or L. Since there is no 'J' in the
					// BLOSUM62 matrix, approximate J = I
					case 'J': case 'j':
					case 'I': case 'i':
						*o = AA::I;
						break;
					case 'L': case 'l':
						*o = AA::L;
						break;
					case 'K': case 'k':
						*o = AA::K;
						break;
					case 'M': case 'm':
						*o = AA::M;
						break;
					case 'F': case 'f':
						*o = AA::F;
						break;
					case 'P': case 'p':
						*o = AA::P;
						break;
					case 'S': case 's':
						*o = AA::S;
						break;
					case 'T': case 't':
						*o = AA::T;
						break;
					case 'W': case 'w':
						*o = AA::W;
						break;
					case 'Y': case 'y':
						*o = AA::Y;
						break;
					case 'V': case 'v':
						*o = AA::V;
						break;
					case 'B': case 'b':
						*o = AA::B;
						break;
					case 'Z': case 'z':
						*o = AA::Z;
						break;
					case 'X': case 'x':
					case 'U': case 'u': // Treat selenocysteine as 'X' (same as BLAST)
						*o = AA::X;
						break;
					default:
						throw __FILE__ ":translate_aa_seq: Illegal base";
						break;
				};
			}
		};
		
		void reverse_complement(std::vector<base_type> &m_seq) const
		{
			#ifdef _DEBUG
			if(!is_na){
				throw __FILE__ ":reverse_complement: Only valid for NA sequences";
			}
			#endif // _DEBUG
			
			size_t len = m_seq.size();

			len = len/2 + (len%2);

			std::vector<base_type>::iterator f = m_seq.begin();
			std::vector<base_type>::reverse_iterator r = m_seq.rbegin();

			for(size_t i = 0;i < len;i++, f++, r++){
				
				const base_type tmp = complement(*f);
				*f = complement(*r);
				*r = tmp;
			}
		};

		inline base_type complement(const base_type &m_base) const
		{
			#ifdef _DEBUG
			if(!is_na){
				throw __FILE__ ":complement: Only valid for NA sequences";
			}
			#endif // _DEBUG
			
			switch(m_base){
				case NA::A:
					return NA::T;
				case NA::C:
					return NA::G;
				case NA::G:
					return NA::C;
				case NA::T:
					return NA::A;
				case NA::M:
					return NA::K;
				case NA::R:
					return NA::Y;
				case NA::S:
					return NA::S;
				case NA::V:
					return NA::B;
				case NA::W:
					return NA::W;
				case NA::Y:
					return NA::R;
				case NA::H:
					return NA::D;
				case NA::K:
					return NA::M;
				case NA::D:
					return NA::H;
				case NA::B:
					return NA::V;
				case NA::N:
					return NA::N;
				case NA::GAP:
					return NA::GAP;
			};

			throw __FILE__ ":complement: Unknown base";
			return NA::GAP; // Keep the compiler happy
		};

		void align_overlap();
		void align_smith_waterman();
		
		SO_Score trace_back();

	public:
		
		SeqOverlap(const AlignmentMode &m_mode, const bool &m_is_na);
			
		~SeqOverlap()
		{
			if(last_row != NULL){
				SIMD_FREE(last_row);
			}
			
			if(curr_row != NULL){
				SIMD_FREE(curr_row);
			}
		};

		void align()
		{
			switch(mode){
				case Overlap:
					align_overlap();
					break;
				case SmithWaterman:
					align_smith_waterman();
					break;
				default:
					throw __FILE__ ":align: Unknown mode";
			};
		}
		
		inline void set_mode(const AlignmentMode &m_mode)
		{
			mode = m_mode;
		};
		
		inline void enable_mask_N_na()
		{
			mask_N_na = true;
		};
		
		inline void disable_mask_N_na()
		{
			mask_N_na = false;
		};
		
		// Pack the supplied query into slot m_index
		inline void pack_query(const unsigned int &m_index, const std::string &m_query)
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":pack_query: m_index >= SO_LEN";
			}
			#endif // _DEBUG
			
			const unsigned int len = m_query.size();

			if(len > MAX_SEQUENCE_LENGTH){
				throw __FILE__ ":pack_query: |m_query| > MAX_SEQUENCE_LENGTH";
			}
			
			if(is_na){
				for(unsigned int i = 0;i < len;i++){
					query[i].v[m_index] = na_to_bits(m_query[i]);
				}
			}
			else{ // aa
				for(unsigned int i = 0;i < len;i++){
					query[i].v[m_index] = aa_to_bits(m_query[i]);
				}
			}

			query_len.v[m_index] = len;
			
			// Recompute the maximum query length
			max_query_len = 0;
			
			for(unsigned int i = 0;i < SO_LEN;i++){			
				max_query_len = std::max(max_query_len, query_len.v[i]);
			}
		};
		
		// Pack the supplied query into *all* slots
		inline void pack_query(const std::string &m_query)
		{
			const unsigned int len = m_query.size();

			if(len > MAX_SEQUENCE_LENGTH){
				throw __FILE__ ":pack_query: |m_query| > MAX_SEQUENCE_LENGTH";
			}
			
			// Update the maximum query length
			max_query_len = len;
			
			if(is_na){
				for(unsigned int i = 0;i < len;i++){
					query[i].sse = SIMD_SET1( na_to_bits(m_query[i]) );
				}
			}
			else{ // aa
				for(unsigned int i = 0;i < len;i++){
					query[i].sse = SIMD_SET1( aa_to_bits(m_query[i]) );
				}
			}
			
			query_len.sse = SIMD_SET1(len);
		};
		
		// Clear all queries in all slots
		inline void clear_query()
		{
			// Reset the maximum query length
			max_query_len = 0;
			
			query_len.sse = SIMD_SET1(0);
		};
		
		// Clear the query in slot m_index
		inline void clear_query(const unsigned int &m_index)
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":clear_query: m_index >= SO_LEN";
			}
			#endif // _DEBUG
			
			query_len.v[m_index] = 0;
			
			// Recompute the maximum query length
			max_query_len = 0;
			
			for(unsigned int i = 0;i < SO_LEN;i++){			
				max_query_len = std::max(max_query_len, query_len.v[i]);
			}
		};
		
		// Pack the supplied target into slot m_index
		inline void pack_target(const unsigned int &m_index, const std::string &m_target)
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":pack_target: m_index >= SO_LEN";
			}
			#endif // _DEBUG
			
			const unsigned int len = m_target.size();
			
			if(len > MAX_SEQUENCE_LENGTH){
				throw __FILE__ ":pack_target: |m_target| > MAX_SEQUENCE_LENGTH";
			}
			
			if(is_na){
				for(unsigned int i = 0;i < len;i++){
					target[i].v[m_index] = na_to_bits(m_target[i]);
				}
			}
			else{ // aa
				for(unsigned int i = 0;i < len;i++){
					target[i].v[m_index] = aa_to_bits(m_target[i]);
				}
			}

			target_len.v[m_index] = len;
			
			// Recompute the maximum target length
			max_target_len = 0;
			
			for(unsigned int i = 0;i < SO_LEN;i++){			
				max_target_len = std::max(max_target_len, target_len.v[i]);
			}
		};
		
		// Pack the supplied target into *all* slots
		inline void pack_target(const std::string &m_target)
		{
			const unsigned int len = m_target.size();

			if(len > MAX_SEQUENCE_LENGTH){
				throw __FILE__ ":pack_target: |m_target| > MAX_SEQUENCE_LENGTH";
			}
			
			// Update the maximum target length
			max_target_len = len;
			
			if(is_na){
				for(unsigned int i = 0;i < len;i++){
					target[i].sse = SIMD_SET1( na_to_bits(m_target[i]) );
				}
			}
			else{ // aa
				for(unsigned int i = 0;i < len;i++){
					target[i].sse = SIMD_SET1( aa_to_bits(m_target[i]) );
				}
			}
			
			target_len.sse = SIMD_SET1(len);
		};
		
		// Clear all targets in all slots
		inline void clear_target()
		{
			
			// Reset the maximum target length
			max_target_len = 0;
			
			target_len.sse = SIMD_SET1(0);
		};
		
		// Clear the target in slot m_index
		inline void clear_target(const unsigned int &m_index)
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":clear_target: m_index >= SO_LEN";
			}
			#endif // _DEBUG
			
			target_len.v[m_index] = 0;
			
			// Recompute the maximum target length
			max_target_len = 0;
			
			for(unsigned int i = 0;i < SO_LEN;i++){
				max_target_len = std::max(max_target_len, target_len.v[i]);
			}
		};
		
		// Return the query and target sequences as std::strings
		std::string query_seq(const unsigned int &m_index) const
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":query_seq: Index out of bounds";
			}
			#endif // _DEBUG
			
			std::string ret;
			
			if(is_na){
				for(SO_Score i = 0;i < query_len.v[m_index];i++){
					ret.push_back( bits_to_na(query[i].v[m_index]) );
				}
			}
			else{
				for(SO_Score i = 0;i < query_len.v[m_index];i++){
					ret.push_back( bits_to_aa(query[i].v[m_index]) );
				}
			}			
			
			return ret;
		};
		
		std::string target_seq(const unsigned int &m_index) const
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":target_seq: Index out of bounds";
			}
			#endif // _DEBUG
			
			std::string ret;
			
			if(is_na){
				for(SO_Score i = 0;i < target_len.v[m_index];i++){
					ret.push_back( bits_to_na(target[i].v[m_index]) );
				}
			}
			else{
				for(SO_Score i = 0;i < target_len.v[m_index];i++){
					ret.push_back( bits_to_aa(target[i].v[m_index]) );
				}
			}		
			
			return ret;
		};
					
		// The coordinates of the first and last aligned base in the query
		inline std::pair<int, int> alignment_range_query(const unsigned int &m_index) const
		{
			
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":alignment_range_query: Index out of bounds";
			}
			#endif // _DEBUG
			
			return std::make_pair(max_elem.M_start_i.v[m_index], stop_i.v[m_index]);
		};
		
		// The coordinates of the first and last aligned base in the query
		inline std::pair<int, int> alignment_range_target(const unsigned int &m_index) const
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":alignment_range_target: Index out of bounds";
			}
			#endif // _DEBUG
			
			return std::make_pair(max_elem.M_start_j.v[m_index], stop_j.v[m_index]);
		};
		
		inline void alignment_range(std::pair<int, int> &m_query_range,
			std::pair<int, int> &m_target_range, const unsigned int &m_index) const
		{
		
			m_query_range = alignment_range_query(m_index);
			m_target_range = alignment_range_target(m_index);
		};
		
		inline SO_Score score(const unsigned int &m_index) const
		{
			#ifdef _DEBUG
			if(m_index >= SO_LEN){
				throw __FILE__ ":score: Index out of bounds";
			}
			#endif // _DEBUG
			
			return max_elem.M.v[m_index];
		};		
};

} // namespace::SO

#endif // __SEQ_OVERLAP
