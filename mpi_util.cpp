#include "mpi_util.h"
#include "motif.h"
#include "genome_difference.h"
#include "motif_options.h"
#include <string.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for std::string
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<string>(const string &m_str)
{
	return sizeof(size_t) + m_str.size();
}

template<>
unsigned char* mpi_pack<string>(unsigned char* m_ptr, const string &m_str)
{
	size_t len = m_str.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	memcpy(m_ptr, m_str.c_str(), len);
	m_ptr += len;
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<string>(unsigned char* m_ptr, string &m_str)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_str.assign( (char*)m_ptr, len );
	m_ptr += len;
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for __uint128_t
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<__uint128_t>(const __uint128_t &m_obj)
{
	return sizeof(__uint128_t);
}

template<>
unsigned char* mpi_pack<__uint128_t>(unsigned char* m_ptr, const __uint128_t &m_obj)
{
	memcpy( m_ptr, &m_obj, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<__uint128_t>(unsigned char* m_ptr, __uint128_t &m_obj)
{
	memcpy( &m_obj, m_ptr, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Motif
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<Motif>(const Motif &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		MOTIF_MEMBERS
	#undef VARIABLE

	return ret;
}

template<>
unsigned char* mpi_pack<Motif>(unsigned char* m_ptr, const Motif &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
    	MOTIF_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

template<>
unsigned char* mpi_unpack<Motif>(unsigned char* m_ptr, Motif &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
    	MOTIF_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for MotifOptions
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<MotifOptions>(const MotifOptions &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		MOTIF_OPTIONS_MEMBERS
	#undef VARIABLE

	return ret;
}

template<>
unsigned char* mpi_pack<MotifOptions>(unsigned char* m_ptr, const MotifOptions &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
    	MOTIF_OPTIONS_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

template<>
unsigned char* mpi_unpack<MotifOptions>(unsigned char* m_ptr, MotifOptions &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
    	MOTIF_OPTIONS_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for GenomeDifference
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<GenomeDifference>(const GenomeDifference &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		GENOME_DIFFERENCE_MEMBERS
	#undef VARIABLE

	return ret;
}

template<>
unsigned char* mpi_pack<GenomeDifference>(unsigned char* m_ptr, const GenomeDifference &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
    	GENOME_DIFFERENCE_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

template<>
unsigned char* mpi_unpack<GenomeDifference>(unsigned char* m_ptr, GenomeDifference &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
    	GENOME_DIFFERENCE_MEMBERS
	#undef VARIABLE

    return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for SequenceBase
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<SequenceBase>(const SequenceBase &m_obj)
{
	return sizeof(m_obj.buffer);
}

template<>
unsigned char* mpi_pack<SequenceBase>(unsigned char* m_ptr, const SequenceBase &m_obj)
{
	m_ptr = mpi_pack(m_ptr, m_obj.buffer);

    return m_ptr;
}

template<>
unsigned char* mpi_unpack<SequenceBase>(unsigned char* m_ptr, SequenceBase &m_obj)
{
	m_ptr = mpi_unpack(m_ptr, m_obj.buffer);

    return m_ptr;
}

