#include "split.h"
#include <sstream>

using namespace std;

// Split input string on white space
vector<string> split(const string &m_line)
{
        vector<string> ret;

        stringstream ssin(m_line);

        string buffer;

        while( ssin >> buffer ){
                ret.push_back(buffer);
        }

        return ret;
}

