#ifndef __MOTIF_OPTIONS
#define __MOTIF_OPTIONS

#include <cstddef> // For size_t
#include <string>
#include <deque>

// There are multiple schemes for dividing the input data into test and training
// sets for cross validation
enum{
	TEST_ON_BRANCH_TIPS,
	TEST_BY_DISTANCE_FROM_ROOT
};

// Forward declaration of mpi helper functions to keep the compiler happy
template <class T> size_t mpi_size(const T &m_obj);
template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);

struct MotifOptions
{
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized.
	#define MOTIF_OPTIONS_MEMBERS \
		VARIABLE(std::string, reference_filename) \
		VARIABLE(std::string, genomes_filename) \
		VARIABLE(std::string, genealogy_filename) \
		VARIABLE(std::string, codon_filename) \
		VARIABLE(std::deque<std::string>, eval_filename) \
		VARIABLE(unsigned int, motif_len) \
		VARIABLE(unsigned int, num_motif) \
		VARIABLE(unsigned int, num_batch) \
		VARIABLE(unsigned int, num_epoch) \
		VARIABLE(unsigned int, num_trial) \
		VARIABLE(unsigned int, seed) \
		VARIABLE(int, cross_validation) \
		VARIABLE(unsigned int, distance_from_root) \
		VARIABLE(bool, use_regions) \
		VARIABLE(bool, quit)

	
	#define VARIABLE(A, B) A B;
		MOTIF_OPTIONS_MEMBERS
	#undef VARIABLE
	
	MotifOptions() {};
	void load(int argc, char* argv[]);

	template<class T> friend size_t mpi_size(const T &m_obj);
	template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
	template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
};

template<> size_t mpi_size(const MotifOptions &m_opt);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const MotifOptions &m_opt);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, MotifOptions &m_opt);

#endif // __MOTIF_OPTIONS