#ifndef __CFG_EPS_REFINE_H__
#define __CFG_EPS_REFINE_H__

#include <refinement/CondEpsGen.h>

namespace covenant {

using namespace std;

enum GeneralizationMethod { GREEDY, MAX_GEN };

std::istream& operator>>(std::istream& in, GeneralizationMethod& gen)
{
    std::string token;
    in >> token;
    if (token == "greedy")
        gen = GREEDY;
    else if (token == "max-gen")
        gen = MAX_GEN;
    else 
      throw error("Invalid generalization method");
    return in;
}

std::ostream& operator<<(std::ostream& o, GeneralizationMethod gen)
{
    if (gen == GREEDY)
      o << "greedy";
    else if (gen == MAX_GEN)
      o << "max-gen";
    else 
      throw error("invalid generalization method");
    return o;
}

template<typename EdgeSym>
class EpsRefine
{
 public:
  
  virtual ~EpsRefine() { }
  
  // Refine a language with respect to a witness, following the chosen
  // generalization strategy.
  // 
  // Return a generalized witness that by construction \not \in
  // language of the CFG.
  virtual DFA<EdgeSym> refine(CondEpsGen<EdgeSym>& gen) = 0;
};

} // end namespace covenant

#endif /*__CFG_REFINE_H__*/

