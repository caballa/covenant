#ifndef __PARSE_UTILS_H__
#define __PARSE_UTILS_H__

#include <Common.h>

namespace covenant{

class StrParse 
{
public:
  StrParse(std::string& _str): idx(0), str(_str) { }

  StrParse(const char* _str): idx(0), str(_str) { }
    
  char peek(void) 
  { 
    return empty() ? EOF : str[idx]; 
  }

  void chomp(void) 
  { 
    if (empty()) throw error("during parsing of a CFG production");
    idx++; 
  }  

  void chomp(char c) 
  { 
    if (empty()) 
      throw error("during parsing of a CFG production");
    if (str[idx] != c)
      throw error("during parsing of a CFG production");
    idx++; 
  }

  char pop(void) 
  { 
    if (empty()) throw error("during parsing of a CFG production");
    return str[idx++]; 
  }

  bool empty(void) 
  { 
    return idx >= str.size(); 
  }

  void reset(const char* _str)
  {
    idx = 0;
    str = _str;
  }

  void reset(std::string& _str)
  {
    idx = 0;
    str = _str;
  }

  //protected:
  unsigned int idx;
  std::string str;
}; /* end class StrParse */

// Parser for an input stream (say, cin)
class StreamParse 
{
public:
  StreamParse(std::istream& _str): str(_str), fresh(false) { }

  char peek(void)     
  {
    if(!fresh)
    {
      if(str.eof())
        return EOF;
      str.get(next);
      fresh = true;
      // std::cout << "Read:" << next << std::endl;
    }
    return next;
  }
  
  void chomp(void) 
  { 
    peek(); 
    fresh = false; 
  }

  void chomp(char c)  
  { 
    peek(); 
    if (next != c)
      throw error("during parsing of a CFG production");
    fresh = false; 
  }

  char pop(void) 
  { 
    peek(); 
    fresh = false; 
    return next; 
  }
  
  bool empty(void) 
  { 
    return (!fresh && str.eof() ); 
  }

  //protected:
  // Input stream
  std::istream& str;
  // Next character
  char next;
  // Has the next character been loaded into 'c'
  bool fresh;

}; // end class StreamParse 

} // end namespace covenant 

#endif /*__PARSE_UTILS_H__*/
