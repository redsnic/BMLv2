#include "StringTokenizer.h"

/* small class to tokenize a string also considering null fields */

/**
 * @brief StringTokenizer::StringTokenizer
 * @param str    a string of any length
 * @param sep    a separator string
 */
StringTokenizer::StringTokenizer(std::string str, std::string sep)
{
    stringList=new std::vector<std::string*>(0);
    position = 0;
    std::size_t begin = 0;
    std::size_t end = 0;

    end = str.find(sep);

    while(end != std::string::npos){     // extract each field, considering null ones
        if(end == begin){
            stringList->push_back(new std::string(""));
            begin++;
            end = str.find(sep, begin);
        }else {
            stringList->push_back(new std::string(str.substr(begin, end-begin)));
            begin = end+1;
            end = str.find(sep, begin);
        }
    }

}

StringTokenizer::~StringTokenizer(){

    for(unsigned int i = 0; i<stringList->size() ; i++){
        delete(stringList->operator [](i));
    }
    delete(stringList);

}

/**
 * @brief StringTokenizer::next
 * @return the value of the next field (an empty string if next is called when end() is true)
 */
std::string StringTokenizer::next(){

    if(position<stringList->size()){
        std::string ret = *(stringList->operator [](position));
        position++;
        return ret;
    } else {
        return "";
    }

}

/**
 * @brief StringTokenizer::end
 * @return true if all fields were returned by next
 */
bool StringTokenizer::end(){
    return position<stringList->size();
}
