#ifndef __SIGMA_STAR_ABSTRACT_H__
#define __SIGMA_STAR_ABSTRACT_H__

///////////////////////////////////////////////////////////////////////
// Convert a CFG into a finite automata by producing Sigma*
///////////////////////////////////////////////////////////////////////

#include <CFG.h>
#include <DFA.h>
#include <abstraction/StronglyRegularToFA.h>

namespace covenant {

template<typename EdgeSym>
class SigmaStarAbstract: public Abstract<EdgeSym> 
{
 public:

  DFA<EdgeSym> do_abstraction(CFG &g, const bool is_regular)
  {
    DFA<EdgeSym> sfa (g.getTermFactory ());
    if (is_regular)
    {        
      LOG("verbose", cout << "The CFG is regular thus no need to abstract.\n");
      StronglyRegularToFA<EdgeSym>  conv(g);
      sfa = conv.make_fa();
    }
    else
    {
      LOG("verbose", cout << "Applied Sigma* abstraction.\n");
      State q0 = sfa.state(0);
      sfa.setStart(q0);
      sfa.accept(q0);
      sfa.eps_transition(q0,q0); 
      for (int i=g.alphstart; i<g.alphstart+g.alphsz; i++)
      {
        sfa.transition(q0,EdgeSym::mkTerm(i), q0);
      }
    }

    sfa.eps_elim();
    DFA<EdgeSym> dfa     = sfa.makeDFA(g.alphstart, g.alphsz);
    DFA<EdgeSym> min_dfa = dfa.minimize(g.alphstart, g.alphsz);
    bool is_not_empty    = min_dfa.eliminateDeadStates();
    if (!is_not_empty)
    {
      throw error("Error during conversion from cfg to reg.");
    }
    
    return min_dfa;
  }
};

} // end namespace covenant

#endif /*__SIGMA_STAR_ABSTRACT_H__*/

