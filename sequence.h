#ifndef __SEQUENCE
#define __SEQUENCE

#include <deque>

template <class T> size_t mpi_size(const T &m_obj);
template<class T> unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
template<class T> unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);

class SequenceBase
{
    private:

        // [4 bits for the codon location][4 bits for the base]
        unsigned char buffer; // <-- Just a single unsigned char!
    
    public:

        // 4 bits per base encoding
        enum {
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
            N = (A | T | C | G)
        };
        
        SequenceBase(const unsigned char m_base = 0)
        {
            // Only set the base (not the codon value)
            buffer = m_base;
        };

        inline void set_ascii_base(char m_base)
        {
            switch(m_base){
                case 'A': case 'a':
                    set_base(A);
                    break;
                case 'C': case 'c':
                    set_base(C);
                    break;
                case 'G': case 'g':
                    set_base(G);
                    break;
                case 'T': case 't':
                    set_base(T);
                    break;
                case 'M': case 'm':
                    set_base(M);
                    break;
                case 'R': case 'r':
                    set_base(R);
                    break;
                case 'S': case 's':
                    set_base(S);
                    break;
                case 'V': case 'v':
                    set_base(V);
                    break;
                case 'W': case 'w':
                    set_base(W);
                    break;
                case 'Y': case 'y':
                    set_base(Y);
                    break;
                case 'H': case 'h':
                    set_base(H);
                    break;
                case 'K': case 'k':
                    set_base(K);
                    break;
                case 'D': case 'd':
                    set_base(D);
                    break;
                case 'B': case 'b':
                    set_base(B);
                    break;
                case 'N': case 'n':
                    set_base(N);
                    break;
                default:
                    throw __FILE__ ":SequenceBase: Invalid base";
            };
        };

        inline void set_base(unsigned char m_base)
        {
            buffer = (buffer & 0xF0)/*codon loc*/ | (m_base & 0x0F)/*base*/;
        };

        inline void set_codon(unsigned char m_codon)
        {
            buffer = (m_codon << 4)/*codon loc*/ | (buffer & 0x0F)/*base*/;
        };

        inline void set_base_and_codon(unsigned char m_base, unsigned char m_codon)
        {
            buffer = (m_codon << 4)/*codon loc*/ | (m_base & 0x0F)/*base*/;
        };

        inline unsigned char codon() const
        {
            // The high-order bits of buffer store the codon location.
            // Since zeros will be shifted in, we do not need to mask the
            // return value.
            return buffer >> 4;
        };

        inline unsigned char base() const
        {
            // The low-order bits of buffer store the binary base
            return buffer & 0x0F;
        };

        inline char ascii_base() const
        {
            // The low-order bits of buffer store the binary base
            switch(buffer & 0x0F){
                case A:
                    return 'A';
                case C:
                    return 'C';
                case G:
                    return 'G';
                case T:
                    return 'T';
                case M:
                    return 'M';
                case R:
                    return 'R';
                case S:
                    return 'S';
                case V:
                    return 'V';
                case W:
                    return 'W';
                case Y:
                    return 'Y';
                case H:
                    return 'H';
                case K:
                    return 'K';
                case D:
                    return 'D';
                case B:
                    return 'B';
                case N:
                    return 'N';
                default:
                    throw __FILE__ ":ascii_base: Invalid base";
            };

            // Should never get here
            return '?';
        };

        template<class T> friend size_t mpi_size(const T &m_obj);
        template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
        template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
};

template<> size_t mpi_size(const SequenceBase &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const SequenceBase &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, SequenceBase &m_obj);

typedef std::deque<SequenceBase> Sequence;

#endif //__SEQUENCE