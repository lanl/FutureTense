#ifndef __MAP_KEYS
#define __MAP_KEYS

#include <deque>
#include <unordered_map>
#include <algorithm>

template<class A, class B>
inline std::deque<A> keys(const std::unordered_multimap<A, B> &m_db)
{
	std::deque<A> ret;

	for(typename std::unordered_multimap<A, B>::const_iterator i = m_db.begin();i != m_db.end();++i){
		ret.push_back(i->first);
	}

	// Make the keys unique
	std::sort( ret.begin(), ret.end() );

	ret.erase( std::unique( ret.begin(), ret.end() ), ret.end() );

	return ret;
};

#endif // __MAP_KEYS