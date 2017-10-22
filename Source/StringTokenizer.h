#ifndef STRINGTOKENIZER_H
#include <stdlib.h>
#include <string>
#include <vector>
#define STRINGTOKENIZER_H


class StringTokenizer
{
public:
    StringTokenizer(std::string, std::string sep);
    std::string next();
    bool end();
    virtual ~StringTokenizer();
private:
    std::vector<std::string*>* stringList;
    unsigned int position;
};

#endif // STRINGTOKENIZER_H
