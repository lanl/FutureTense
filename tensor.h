#ifndef __TENSOR
#define __TENSOR

#include <vector>

template<class T>
class Tensor : public std::vector<T>
{
    public:
        typedef unsigned int size_type;

    private:

        // Dimensionality of tensor = length.size()
        // Number of elements in i^th dimension = length[i]
        std::vector<size_type> length;

        inline size_t product(const typename std::vector<size_type>::const_iterator &m_begin,
            const typename std::vector<size_type>::const_iterator &m_end) const
        {
            size_t ret = 1;

            for(typename std::vector<size_type>::const_iterator i = m_begin;i != m_end;++i){
                ret *= *i;
            }

            return ret;
        };

    public:

        Tensor()
        {

        };

        Tensor(const std::vector<size_type> &m_dim) :
            length(m_dim)
        {
            std::vector<T>::resize( product(length.begin(), length.end()) );
        };

        Tensor( const size_type &m_a)
        {
            length.resize(1);

            length[0] = m_a;

            std::vector<T>::resize( product(length.begin(), length.end()) );
        };

        Tensor(const size_type &m_a, const size_type &m_b)
        {
            length.resize(2);

            length[0] = m_a;
            length[1] = m_b;

            std::vector<T>::resize( product(length.begin(), length.end()) );
        };

        Tensor(const size_type &m_a, const size_type &m_b, const size_type &m_c)
        {
            length.resize(3);

            length[0] = m_a;
            length[1] = m_b;
            length[2] = m_c;

            std::vector<T>::resize( product(length.begin(), length.end()) );
        };

        Tensor(const size_type &m_a, const size_type &m_b, const size_type &m_c, const size_type &m_d)
        {
            length.resize(4);

            length[0] = m_a;
            length[1] = m_b;
            length[2] = m_c;
            length[3] = m_d;

            std::vector<T>::resize( product(length.begin(), length.end()) );
        };

        Tensor(const size_type &m_a, const size_type &m_b, const size_type &m_c, const size_type &m_d, const size_type &m_e)
        {
            length.resize(5);

            length[0] = m_a;
            length[1] = m_b;
            length[2] = m_c;
            length[3] = m_d;
            length[4] = m_e;

            std::vector<T>::resize( product(length.begin(), length.end()) );
        };

        // Allow setting all tensor elements from a single scalar
        inline Tensor<T>& operator=(const T& m_rhs)
        {
            for(typename std::vector<T>::iterator i = std::vector<T>::begin();i != std::vector<T>::end();++i){
                *i = m_rhs;
            }

            return *this;
        };

        // Multiply all tensor elements by a scalar
        inline Tensor<T>& operator*=(const float& m_rhs)
        {
            for(typename std::vector<T>::iterator i = std::vector<T>::begin();i != std::vector<T>::end();++i){
                *i *= m_rhs;
            }

            return *this;
        };

         // Divide all tensor elements by a scalar
        inline Tensor<T>& operator/=(const float& m_rhs)
        {
            for(typename std::vector<T>::iterator i = std::vector<T>::begin();i != std::vector<T>::end();++i){
                *i /= m_rhs;
            }

            return *this;
        };

        // Increment all tensor elements by the tensor elements of an equally sized tensor
        inline Tensor<T>& operator+=(const Tensor<T>& m_rhs)
        {
            #ifdef EBUG
            if( length != m_rhs.length ){
                throw __FILE__ ":Tensor+=: Cannot add tensors with unequal sizes";
            }
            #endif // EBUG

            typename std::vector<T>::const_iterator j = m_rhs.begin();

            for(typename std::vector<T>::iterator i = std::vector<T>::begin();i != std::vector<T>::end();++i, ++j){
                *i += *j;
            }

            return *this;
        };

        inline T& operator()(const size_type &m_a)
        {
            #ifdef EBUG
            if(length.size() != 1){
                throw __FILE__ ":Tensor(1): Index has incorrect dimensionality";
            }

            if(m_a >= length[0]){
                throw __FILE__ ":Tensor(1): Index is out of bounds";
            }
            #endif // EBUG

            return (*this)[m_a];
        };

        inline const T& operator()(const size_type &m_a) const
        {
            #ifdef EBUG
            if(length.size() != 1){
                throw __FILE__ ":Tensor(1): Index has incorrect dimensionality";
            }

            if(m_a >= length[0]){
                throw __FILE__ ":Tensor(1): Index is out of bounds";
            }
            #endif // EBUG

            return (*this)[m_a];
        };

        inline T& operator()(const size_type &m_a, const size_type &m_b)
        {
            #ifdef EBUG
            if(length.size() != 2){
                throw __FILE__ ":Tensor(2): Index has incorrect dimensionality";
            }

            if( (m_a >= length[0]) || (m_b >= length[1]) ){
                throw __FILE__ ":Tensor(2): Index is out of bounds";
            }
            #endif // EBUG

            const size_type master_index = m_a*length[1] + m_b;

            #ifdef EBUG
            if( master_index >= std::vector<T>::size() ){
                throw __FILE__ ":Tensor(2): Master index is out of bounds";
            }
            #endif // EBUG

            return (*this)[master_index];
        };

        inline const T& operator()(const size_type &m_a, const size_type &m_b) const
        {
            #ifdef EBUG
            if(length.size() != 2){
                throw __FILE__ ":Tensor(2): Index has incorrect dimensionality";
            }

            if( (m_a >= length[0]) || (m_b >= length[1]) ){
                throw __FILE__ ":Tensor(2): Index is out of bounds";
            }
            #endif // EBUG

            const size_type master_index = m_a*length[1] + m_b;

            #ifdef EBUG
            if( master_index >= std::vector<T>::size() ){
                throw __FILE__ ":Tensor(2): Master index is out of bounds";
            }
            #endif // EBUG

            return (*this)[master_index];
        };

        inline T& operator()(const size_type &m_a, const size_type &m_b, const size_type &m_c)
        {
            #ifdef EBUG
            if(length.size() != 3){
                throw __FILE__ ":Tensor(3): Index has incorrect dimensionality";
            }

            if( (m_a >= length[0]) || (m_b >= length[1]) || (m_c >= length[2])){
                throw __FILE__ ":Tensor(3): Index is out of bounds";
            }
            #endif // EBUG

            const size_type master_index = (m_a*length[1] + m_b)*length[2] + m_c;

            #ifdef EBUG
            if( master_index >= std::vector<T>::size() ){
                throw __FILE__ ":Tensor(3): Master index is out of bounds";
            }
            #endif // EBUG

            return (*this)[master_index];
        };

        inline const T& operator()(const size_type &m_a, const size_type &m_b, const size_type &m_c) const
        {
            #ifdef EBUG
            if(length.size() != 3){
                throw __FILE__ ":Tensor(3): Index has incorrect dimensionality";
            }

            if( (m_a >= length[0]) || (m_b >= length[1]) || (m_c >= length[2])){
                throw __FILE__ ":Tensor(3): Index is out of bounds";
            }

            #endif // EBUG

            const size_type master_index = (m_a*length[1] + m_b)*length[2] + m_c;

            #ifdef EBUG
            if( master_index >= std::vector<T>::size() ){
                throw __FILE__ ":Tensor(3): Master index is out of bounds";
            }
            #endif // EBUG

            return (*this)[master_index];
        };

        inline T& operator()(const size_type &m_a, const size_type &m_b, const size_type &m_c, const size_type &m_d)
        {
            #ifdef EBUG
            if(length.size() != 4){
                throw __FILE__ ":Tensor(4): Index has incorrect dimensionality";
            }

            if( (m_a >= length[0]) || (m_b >= length[1]) || (m_c >= length[2]) || (m_d >= length[3]) ){
                throw __FILE__ ":Tensor(4): Index is out of bounds";
            }
            #endif // EBUG

            const size_type master_index = ( (m_a*length[1] + m_b)*length[2] + m_c)*length[3] + m_d;

            #ifdef EBUG
            if( master_index >= std::vector<T>::size() ){
                throw __FILE__ ":Tensor(4): Master index is out of bounds";
            }
            #endif // EBUG

            return (*this)[master_index];
        };

        inline const T& operator()(const size_type &m_a, const size_type &m_b, const size_type &m_c, const size_type &m_d) const
        {
            #ifdef EBUG
            if(length.size() != 4){
                throw __FILE__ ":Tensor(4): Index has incorrect dimensionality";
            }

            if( (m_a >= length[0]) || (m_b >= length[1]) || (m_c >= length[2]) || (m_d >= length[3]) ){
                throw __FILE__ ":Tensor(4): Index is out of bounds";
            }

            #endif // EBUG

            const size_type master_index = ( (m_a*length[1] + m_b)*length[2] + m_c)*length[3] + m_d;

            #ifdef EBUG
            if( master_index >= std::vector<T>::size() ){
                throw __FILE__ ":Tensor(4): Master index is out of bounds";
            }
            #endif // EBUG

            return (*this)[master_index];
        };

        inline T& operator()(const size_type &m_a, const size_type &m_b, const size_type &m_c, const size_type &m_d, 
            const size_type &m_e)
        {
            #ifdef EBUG
            if(length.size() != 5){
                throw __FILE__ ":Tensor(5): Index has incorrect dimensionality";
            }

            if( (m_a >= length[0]) || (m_b >= length[1]) || (m_c >= length[2]) || (m_d >= length[3]) || (m_d >= length[4]) ){
                throw __FILE__ ":Tensor(5): Index is out of bounds";
            }
            #endif // EBUG

            const size_type master_index = ( ( (m_a*length[1] + m_b)*length[2] + m_c)*length[3] + m_d)*length[4] + m_e;

            #ifdef EBUG
            if( master_index >= std::vector<T>::size() ){
                throw __FILE__ ":Tensor(5): Master index is out of bounds";
            }
            #endif // EBUG

            return (*this)[master_index];
        };

        inline const T& operator()(const size_type &m_a, const size_type &m_b, const size_type &m_c, const size_type &m_d, 
            const size_type &m_e) const
        {
            #ifdef EBUG
            if(length.size() != 5){
                throw __FILE__ ":Tensor(5): Index has incorrect dimensionality";
            }

            if( (m_a >= length[0]) || (m_b >= length[1]) || (m_c >= length[2]) || (m_d >= length[3]) || (m_d >= length[4]) ){
                throw __FILE__ ":Tensor(5): Index is out of bounds";
            }
            #endif // EBUG

            const size_type master_index = ( ( (m_a*length[1] + m_b)*length[2] + m_c)*length[3] + m_d)*length[4] + m_e;

            #ifdef EBUG
            if( master_index >= std::vector<T>::size() ){
                throw __FILE__ ":Tensor(5): Master index is out of bounds";
            }
            #endif // EBUG

            return (*this)[master_index];
        };

        inline size_t dim() const
        {
            return length.size();
        };

        inline size_t size(const size_type &m_dim) const
        {
            #ifdef EBUG
            if( m_dim >= length.size() ){
                throw __FILE__ ":Tensor::size(): Dimension out of bounds";
            }
            #endif // EBUG

            return length[m_dim];
        };

        inline const std::vector<size_type>& dim_size() const
        {
            return length;
        };

        inline void clear()
        {
            length.clear();

            std::vector<size_type>::clear();
        };

        void normalize()
        {
            T total = 0.0;

            for(typename std::vector<T>::const_iterator i = std::vector<T>::begin();i != std::vector<T>::end();++i){
                total += *i;
            }

            for(typename std::vector<T>::iterator i = std::vector<T>::begin();i != std::vector<T>::end();++i){
                *i /= total;
            }
        };
};

#endif // __TENSOR
