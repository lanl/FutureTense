#ifndef __GENOME_DIFFERENCE
#define __GENOME_DIFFERENCE

#include <deque>
#include <sstream>

// We need to protect any commas that appear in template variables
// if we are going to combine them with X Macros
#define SINGLE_ARG(...) __VA_ARGS__

template <class T> size_t mpi_size(const T &m_obj);
template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);

class GenomeDifference
{
	private:

		// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
		// ensure that structure variable are correctly serialized.
		#define GENOME_DIFFERENCE_MEMBERS \
			VARIABLE(SINGLE_ARG(std::deque< std::pair<unsigned int, char> >), substitution) \
			VARIABLE(SINGLE_ARG(std::deque< std::pair<unsigned int, std::string> >), insertion) \
			VARIABLE(std::deque<unsigned int>, deletion)

		#define VARIABLE(A, B) A B;
			GENOME_DIFFERENCE_MEMBERS
		#undef VARIABLE

		std::string format_insertion(const std::string &m_str) const
		{
			std::string ret = m_str;

			for(std::string::iterator i = ret.begin();i != ret.end();++i){
				*i = tolower(*i);
			}

			return ret;
		};

	public:

		GenomeDifference()
		{
			// Do nothing
		};

		GenomeDifference(const std::string &m_buffer)
		{
			*this = m_buffer;
		}

		GenomeDifference& operator=(const std::string &m_rhs)
		{
			// An empty string means that this genome is identical to the reference
			if( m_rhs.empty() ){
				return *this;
			}

			if( (m_rhs.front() < '0') || (m_rhs.front() > '9') ){
				throw __FILE__ ":GenomeDifference::operator=: Malformed input string";
			}

			size_t number = 0;
			std::string insertion_buffer;

			for(std::string::const_iterator i = m_rhs.begin();i != m_rhs.end();++i){

				// Since we are using '-' for the gap symbol, we can't use isdigit()!
				if( (*i >= '0') && (*i <= '9') ){

					if( !insertion_buffer.empty() ){

						insertion.push_back( std::make_pair(number, insertion_buffer) );

						number = 0;
						insertion_buffer.clear();
					}

					number = 10*number + (*i - '0');
				}
				else{ // Not a digit

					// Inserted sequences are lower case
					if( islower(*i) ){
						insertion_buffer.push_back(*i);
					}
					else{

						if(*i == '-'){
							deletion.push_back(number);
						}
						else{
							substitution.push_back( std::make_pair(number, *i) );
						}

						number = 0;
					}
				}
			}

			if( !insertion_buffer.empty() ){
				insertion.push_back( std::make_pair(number, insertion_buffer) );
			}

			return *this;
		};

		inline unsigned int edit_distance_to_reference() const
		{
			unsigned int ret = substitution.size() + deletion.size();

			for(std::deque< std::pair<unsigned int, std::string> >::const_iterator i = insertion.begin();i != insertion.end();++i){
				ret += i->second.size();
			}

			return ret;
		};

		// Compute the edit distance between two genomes.
		// Note that (a - b) = (b - a)
		unsigned int distance(const GenomeDifference &m_rhs) const
		{
			unsigned int ret = 0;

			/////////////////////////////////////////////////////////
			// Substitutions
			/////////////////////////////////////////////////////////
			std::deque< std::pair<unsigned int, char> >::const_iterator sub_a = substitution.begin();
			std::deque< std::pair<unsigned int, char> >::const_iterator sub_b = m_rhs.substitution.begin();

			while( ( sub_a != substitution.end() ) || (  sub_b != m_rhs.substitution.end()  ) ){

				if( sub_a == substitution.end() ){

					++sub_b;
					++ret;
					continue;
				}

				if( sub_b == m_rhs.substitution.end() ){

					++sub_a;
					++ret;
					continue;
				}

				if(sub_a->first < sub_b->first){

					++sub_a;
					++ret;
					continue;
				}

				if(sub_a->first > sub_b->first){

					++sub_b;
					++ret;
					continue;
				}

				// sub_a->first == sub_b->first
				ret += (sub_a->second == sub_b->second) ? 0 : 1;
				++sub_a;
				++sub_b;
			}

			// DEBUG
			//cout << "ret(sub) = " << ret << endl;

			/////////////////////////////////////////////////////////
			// Insertions
			/////////////////////////////////////////////////////////
			std::deque< std::pair<unsigned int, std::string> >::const_iterator ins_a = insertion.begin();
			std::deque< std::pair<unsigned int, std::string> >::const_iterator ins_b = m_rhs.insertion.begin();

			while( ( ins_a != insertion.end() ) || (  ins_b != m_rhs.insertion.end()  ) ){

				if( ins_a == insertion.end() ){

					ret += ins_b->second.size();
					++ins_b;
					continue;
				}

				if( ins_b == m_rhs.insertion.end() ){

					ret += ins_a->second.size();
					++ins_a;
					continue;
				}

				if(ins_a->first < ins_b->first){

					ret += ins_a->second.size();
					++ins_a;
					continue;
				}

				if(ins_a->first > ins_b->first){

					ret += ins_b->second.size();
					++ins_b;
					continue;
				}

				// For now, consider non-identical inserted sequences as being completely distinct
				// ins_a->first == ins_b->first
				ret += (ins_a->second == ins_b->second) ? 0 : 
					ins_a->second.size() + ins_b->second.size();

				++ins_a;
				++ins_b;
			}

			// DEBUG
			//cout << "ret(ins) = " << ret << endl;

			/////////////////////////////////////////////////////////
			// deletions
			/////////////////////////////////////////////////////////
			std::deque<unsigned int>::const_iterator del_a = deletion.begin();
			std::deque<unsigned int>::const_iterator del_b = m_rhs.deletion.begin();

			while( ( del_a != deletion.end() ) || (  del_b != m_rhs.deletion.end()  ) ){

				if( del_a == deletion.end() ){

					++del_b;
					++ret;
					continue;
				}

				if( del_b == m_rhs.deletion.end() ){

					++del_a;
					++ret;
					continue;
				}

				if(*del_a < *del_b){

					++del_a;
					++ret;
					continue;
				}

				if(*del_a > *del_b){

					++del_b;
					++ret;
					continue;
				}

				ret += (*del_a == *del_b) ? 0 : 1;
				++del_a;
				++del_b;
			}

			return ret;
		};

		std::string str() const
		{
			std::stringstream ssout;

			std::deque< std::pair<unsigned int, char> >::const_iterator s = substitution.begin();
			std::deque< std::pair<unsigned int, std::string> >::const_iterator i = insertion.begin();
			std::deque<unsigned int>::const_iterator d = deletion.begin();

			const unsigned int max_loc = 0xFFFFFFFF;

			while( ( s != substitution.end() ) || ( i != insertion.end() ) || ( d != deletion.end() ) ){

				const unsigned int s_loc = ( s == substitution.end() ) ? max_loc : s->first;
				const unsigned int i_loc = ( i == insertion.end() ) ? max_loc : i->first;
				const unsigned int d_loc = ( d == deletion.end() ) ? max_loc : *d;

				if( (s_loc <= i_loc) && (s_loc <= d_loc) ){

					ssout << s->first << s->second;
					++s;
					continue;
				}

				if( (i_loc <= s_loc) && (i_loc <= d_loc) ){

					ssout << i->first << format_insertion(i->second);
					++i;
					continue;
				}

				if( (d_loc <= s_loc) && (d_loc <= i_loc) ){

					ssout << *d << '-';
					++d;
					continue;
				}
			}

			return ssout.str();
		};

		const std::deque< std::pair<unsigned int, char> >& get_substitution() const
		{
			return substitution;
		};

		const std::deque< std::pair<unsigned int, std::string> >& get_insertion() const
		{
			return insertion;
		};

		const std::deque<unsigned int>& get_deletion() const
		{
			return deletion;
		};

		bool is_ATGC() const
		{
			// Assuming that the reference genome only contains ATGC, check the substitutions and insertions
			// for any non-ATGC bases

			for(std::deque< std::pair<unsigned int, char> >::const_iterator i = substitution.begin();i != substitution.end();++i){
				switch(i->second){
					case 'A': case 'a':
					case 'T': case 't':
					case 'G': case 'g':
					case 'C': case 'c':
						break;
					default:
						return false;
				};
			}

			for(std::deque< std::pair<unsigned int, std::string> >::const_iterator i = insertion.begin();i != insertion.end();++i){

				for(std::string::const_iterator j = i->second.begin();j != i->second.end();++j){
					switch(*j){
						case 'A': case 'a':
						case 'T': case 't':
						case 'G': case 'g':
						case 'C': case 'c':
							break;
						default:
							return false;
					};
				}
			}

			return true;
		};

		// Recover the original sequence string
		std::string str(const std::string &m_seq) const
		{
			std::deque<char> ret;

			const unsigned int len = m_seq.size();

			std::deque< std::pair<unsigned int, char> >::const_iterator sub_iter = substitution.begin();
			std::deque< std::pair<unsigned int, std::string> >::const_iterator ins_iter = insertion.begin();
			std::deque<unsigned int>::const_iterator del_iter = deletion.begin();

			for(unsigned int i = 0;i < len;++i){
				
				if( ( sub_iter != substitution.end() ) && (i == sub_iter->first) ){

					ret.push_back( toupper(sub_iter->second) );
					++sub_iter;
					continue;
				}

				if( ( ins_iter != insertion.end() ) && (i == ins_iter->first) ){

					for(std::string::const_iterator j = ins_iter->second.begin();j != ins_iter->second.end();++j){
						ret.push_back( toupper(*j) );
					}
					
					++ins_iter;
					continue;
				}

				if( ( del_iter != deletion.end() ) && (i == *del_iter) ){

					// Don't add any base to the output sequence!
					++del_iter;
					continue;
				}

				ret.push_back( toupper(m_seq[i]) );
			}

			return std::string( ret.begin(), ret.end() );
		};

		template<class T> friend size_t mpi_size(const T &m_obj);
        template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
        template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
};

template<> size_t mpi_size(const GenomeDifference &m_diff);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const GenomeDifference &m_diff);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, GenomeDifference &m_diff);

#endif // __GENOME_DIFFERENCE