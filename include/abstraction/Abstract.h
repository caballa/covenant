#ifndef __CFG_ABSTRACT_H__
#define __CFG_ABSTRACT_H__

#include <iostream>
#include <CFG.h>
#include <DFA.h>

namespace covenant {

enum AbstractMethod { SIGMA_STAR, CYCLE_BREAKING};

std::istream& operator>>(std::istream& in, AbstractMethod& abs)
{
    std::string token;
    in >> token;
    if (token == "sigma-star")
        abs = SIGMA_STAR;
    else if (token == "cycle-breaking")
        abs = CYCLE_BREAKING;
    else throw error("invalid abstraction method");
    return in;
}

std::ostream& operator<<(std::ostream& o, AbstractMethod abs)
{
    if (abs == SIGMA_STAR)
      o << "sigma-star";
    else if (abs == CYCLE_BREAKING)
      o << "cycle-breaking";
    else throw error("invalid abstraction method");
    return o;
}


template<typename EdgeSym>
class Abstract{

 public:

  virtual ~Abstract(){}
  
  // Given a CFG g generates a finite automata a such that L(g) is
  // contained in L(a).
  virtual DFA<EdgeSym> 
  do_abstraction(const CFG &g, const bool is_regular) = 0;

};

} // end namespace covenant
#endif /*__CFG_ABSTRACT_H__*/
