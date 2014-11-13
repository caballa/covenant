#ifndef __CFG_GREEDY_EPS_REFINE_H__
#define __CFG_GREEDY_EPS_REFINE_H__

#include <refinement/EpsRefine.h>

namespace covenant {

  using namespace std;

  template<typename EdgeSym>
  class GreedyEpsRefine : public EpsRefine<EdgeSym> 
  {
    
   public:

    typedef DFA<EdgeSym> dfa_t;
    typedef CondEpsGen<EdgeSym> cond_eps_gen_t;


    // Refine a language with respect to a witness, following the
    // chosen generalization strategy.
    // 
    // Return a generalized witness that by construction \not \in
    // language of the CFG.
    dfa_t refine(cond_eps_gen_t& gen) 
    {
      if(gen.word.size() == 0)
      {
        // Just give an empty automaton.
        dfa_t dfa (gen.getTermFactory ());
        State q0 = dfa.state(0);
        dfa.setStart(q0);
        dfa.accept(q0);
        return dfa;
      }
      // Try greedily adding backwards transitions first. 
      for(int di = 0; di < gen.nStates(); di++)
      {
        for(int qi = 0; qi < gen.nStates() - di; qi++)
        {
          int qj = qi + di;
          // Push a choice point
          gen.push();
          if(!gen.addTrans(qj, gen.word[qj], qi))
            gen.pop();
        }
      }
      
      // Now add forward-epsilon transitions
      for(int di = 1; di < gen.nStates(); di++)
      {
        for(int qi = 0; qi < gen.nStates() - di; qi++)
        {
          int qj = qi + di;
          // Push a choice point
          gen.push();
          if(!gen.addEpsTrans(qi, qj))
            gen.pop();
        }
      }

      return gen.emit_dfa();
    }
  };

} // end namespace covenant

#endif /*__CFG_REFINE_H__*/

