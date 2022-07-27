#ifndef __MOTIF
#define __MOTIF

#include <vector>
#include <ostream>
#include <algorithm>
#include <math.h>
#include "word.h"
#include "sequence.h"

#define		MOTIF_VERSION		"0.98 -- Improved objective function"
//	- Improved the objective function to include all of the alternate mutational variants at at given position
//	  when computing the AUROC. The previous version only included the highest scoring variant at a given position.
//	  This lead to all of the variant scores for a parent kmer being set to the same value.
//
//#define		MOTIF_VERSION		"0.97 -- Improved null model"
//	- Improved the null model to select the most frequent substitution for each parent base, codon and region.
//
//#define		MOTIF_VERSION		"0.96 -- Explore optimization strategies"
//	- Changed the rescale_motif() function to make the smallest value 0 and the largest value 1
//
//#define		MOTIF_VERSION		"0.95 -- Better optimization"
//	- Changed the default number of trials to 1. Optimization appears to be more reliable when starting from
//	  all motif elements initialized to zero instead of randomly sampled from a uniform distribution. Only randomly
//	  sample starting motif values when performing more than one optimization trial.
//		- Note that initializing all motif element values to the same value (for example, 0) yields an equal score for
//		  all mutations and therfore a AUROC value of 0.5. I tried randomly assigning starting motif values from a 
//		  uniform distribution [0,1] which would frequently have a starting AUROC value < 0.5.
//	- Fixed the bug of failing to synchronize the test and training lists amoung all MPI ranks. This effected
// 	  the final calculation of the test set AUROC!
//	- Changed the default number of batches to 1.
//
//#define		MOTIF_VERSION		"0.9 -- Code optimization"
//	- Minor code tweaks to speed up implementation.
//	- Cleaned up the implementations of Motif and WordIndex to use less memory
//
//#define		MOTIF_VERSION		"0.8 -- Numeric limits"
//	- Fixed use of std::numeric_limits
//
//#define		MOTIF_VERSION		"0.7 -- Codon locations"
//	- Added location within a codon as a conditional variable
//
//#define		MOTIF_VERSION		"0.6 -- branch length based cross validation"
//#define		MOTIF_VERSION		"0.5 -- Fixed ensemble average bug"
//#define		MOTIF_VERSION		"0.4 -- optional region params"
//#define		MOTIF_VERSION		"0.3 -- Ensemble learning"
//#define		MOTIF_VERSION		"0.2 -- region params"
//#define		MOTIF_VERSION		"0.1 -- In the beginning ..."

// The use of an unsigned short int for storing the word index in WordIndex limits us to
// a maximum motif length of 8, since 4^8 = 2^16
//#define		MAX_MOTIF_LEN		8
// The use of an unsigned int for storing the word index in WordIndex limits us to
// a maximum motif length of 16, since 4^16 = 2^32
#define		MAX_MOTIF_LEN		16

template <class T> size_t mpi_size(const T &m_obj);
template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);

class Motif
{
	private:

		// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
		// ensure that structure variable are correctly serialized.
		#define MOTIF_MEMBERS \
			VARIABLE(std::vector<float>, weight) \
			VARIABLE(unsigned int, motif_len)

		#define VARIABLE(A, B) A B;
			MOTIF_MEMBERS
		#undef VARIABLE

	public:

		typedef std::vector<float>::iterator iterator;
		typedef std::vector<float>::const_iterator const_iterator;

		Motif(const unsigned int &m_len = 0, const float &m_value = 0.0f) : motif_len(m_len)
		{
			weight.resize(NUM_BASE*motif_len, m_value);
		};

		inline iterator begin()
		{
			return weight.begin();
		};

		inline const_iterator begin() const
		{
			return weight.begin();
		};

		inline iterator end()
		{
			return weight.end();
		};

		inline const_iterator end() const
		{
			return weight.end();
		};

		inline unsigned int size() const
		{
			return motif_len;
		};

		inline float operator()(Word m_seq /* copy */) const
		{
			#ifdef EBUG
			if(motif_len == 0){
				throw __FILE__ ":Motif(): Can't score a zero-length motif";
			}
			#endif // EBUG

			float ret = 0.0;

			// The fragment score is evaluated in reverse order (from 3' -> 5')
			unsigned int index = motif_len - 1;

			for(unsigned int i = 0;i < motif_len;++i, --index, m_seq >>= 2){

				// Using two bits per base, mask each base with b11 = 3
				ret += (*this)(m_seq & 3 /*base*/, index /*position*/); // <-- reverse order in position
			}

			return ret;
		};

		inline float& operator()(const unsigned int &m_base /*row*/, const unsigned int &m_index /*column*/)
		{
			#ifdef EBUG
			if( m_base >= NUM_BASE ){
				throw __FILE__ ":Motif(): base out of bounds";
			}

			if(m_index >= motif_len){
				throw __FILE__ ":Motif(): Index out of bounds";
			}
			#endif // EBUG
			
			return weight[m_index*NUM_BASE + m_base];
		};

		inline float operator()(const unsigned int &m_base /*row*/, const unsigned int &m_index /*column*/) const
		{
			#ifdef EBUG
			if( m_base >= NUM_BASE ){
				throw __FILE__ ":Motif(): base out of bounds";
			}

			if(m_index >= motif_len){
				throw __FILE__ ":Motif(): Index out of bounds";
			}
			#endif // EBUG
			
			return weight[m_index*NUM_BASE + m_base];
		};

		// Initialize all elements to the same value
		inline Motif& operator=(const float &m_value)
		{
			for(std::vector<float>::iterator i = weight.begin();i != weight.end();++i){
				*i = m_value;
			}

			return *this;
		};

		// Scale the motif by a constant
		inline Motif operator*(const float &m_scale) const
		{
			Motif ret(*this);

			ret *= m_scale;

			return ret;
		};

		inline Motif& operator*=(const float &m_value)
		{
			for(std::vector<float>::iterator i = weight.begin();i != weight.end();++i){
				*i *= m_value;
			}

			return *this;
		};

		inline Motif& operator/=(const float &m_value)
		{
			for(std::vector<float>::iterator i = weight.begin();i != weight.end();++i){
				*i /= m_value;
			}

			return *this;
		};

		// Add motifs (useful for gradient-based optimization)
		inline Motif& operator+=(const Motif &m_rhs)
		{
			if(motif_len == 0){

				*this = m_rhs;
				return *this;
			}

			#ifdef EBUG
			if(motif_len != m_rhs.motif_len){
				throw __FILE__ ":Motif::operator+=(): |lhs| != |rhs|";
			}
			#endif // EBUG

			const unsigned int total_len = weight.size();

			for(unsigned int i = 0;i < total_len;++i){
				weight[i] += m_rhs.weight[i];
			}

			return *this;
		};

		inline Motif& operator-=(const Motif &m_rhs)
		{
			if(motif_len == 0){
				*this = Motif(m_rhs.motif_len);
			}

			#ifdef EBUG
			if(motif_len != m_rhs.motif_len){
				throw __FILE__ ":Motif::operator-=(): |lhs| != |rhs|";
			}
			#endif // EBUG

			const unsigned int total_len = weight.size();

			for(unsigned int i = 0;i < total_len;++i){
				weight[i] -= m_rhs.weight[i];
			}

			return *this;
		};

		// Return the largest element in the PSWM
		inline float max() const
		{
			return *std::max_element( weight.begin(), weight.end() );
		};

		// Return the smallest element in the PSWM
		inline float min() const
		{
			return *std::min_element( weight.begin(), weight.end() );
		};

		// Return the element with the largest absolute value
		inline float max_abs() const
		{
			float ret = 0.0;

			for(std::vector<float>::const_iterator i = weight.begin();i != weight.end();++i){
				if (fabs(ret) < fabs(*i)){
					ret = *i;
				}
			}

			return ret;
		};

		inline float sum_of_squares() const
		{
			float ret = 0.0;

			for(std::vector<float>::const_iterator i = weight.begin();i != weight.end();++i){
				ret += (*i) * (*i);
			}

			return ret;
		};

		inline float distance(const Motif &m_rhs) const
		{
			#ifdef EBUG
			if(motif_len != m_rhs.motif_len){
				throw __FILE__ ":Motif::distance: size mismatch";
			}
			#endif // EBUG

			float ret = 0.0;

			const unsigned int total_len = weight.size();

			for(unsigned int i = 0;i < total_len;++i){

				const float diff = weight[i] - m_rhs.weight[i];

				ret += diff*diff;
			}

			return sqrt( fabs(ret) );
		};

		template<class T> friend size_t mpi_size(const T &m_obj);
        template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
        template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
};

template<> size_t mpi_size(const Motif &m_motif);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const Motif &m_motif);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, Motif &m_motif);

std::ostream& operator<<(std::ostream& m_out, const Motif& m_motif);

#endif // __MOTIF