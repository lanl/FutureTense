#ifndef __RANDOM_SHUFFLE
#define __RANDOM_SHUFFLE

// Random number generator
#include <gsl/gsl_rng.h>

// A random_shuffle-like function that uses the GSL
template <class T>
void randomize(const T &m_begin, const T &m_end, const gsl_rng *m_ptr_rand)
{
	const size_t len = m_end - m_begin;
	
	for(size_t i = 0;i < len;++i){
	
		const size_t index = gsl_rng_uniform_int(m_ptr_rand, len);
		
		std::swap( *(m_begin + i), *(m_begin + index) );
	}
}

#endif // __RANDOM_SHUFFLE
