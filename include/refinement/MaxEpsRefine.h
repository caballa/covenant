#ifndef __CFG_MAX_EPS_REFINE_H__
#define __CFG_MAX_EPS_REFINE_H__

#include <refinement/EpsRefine.h>

namespace covenant {

  using namespace std;

  template<typename EdgeSym>
  class MaxEpsRefine : public EpsRefine<EdgeSym> {

    int alphstart;
    int alphsz;

  public:

    typedef DFA<EdgeSym> dfa_t;
    typedef CondEpsGen<EdgeSym> cond_eps_gen_t;

    MaxEpsRefine(int _alphstart, int _alphsz)
      : alphstart(_alphstart), alphsz(_alphsz){ } 

    ~MaxEpsRefine() { }

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
        dfa_t dfa;
        State q0 = dfa.state(0);
        dfa.setStart(q0);
        dfa.accept(q0);
        return dfa;
      }
      return maximum_refinement(gen, 0, 0);
    }

    dfa_t maximum_refinement(cond_eps_gen_t& gen, int qi, int qj)
    {
      if(qi >= gen.nStates())
        return gen.emit_dfa();
      assert(qj < gen.nStates());

      int qj_next = (qj+1)%gen.nStates();
      int qi_next = (qj_next == 0) ? qi+1 : qi;
      
      dfa_t ret_dfa(maximum_refinement(gen, qi_next, qj_next));
      ret_dfa.eps_elim();

      gen.push();

      if(qi < qj)
      {
        // Forwards; an epsilon transition
        if(gen.addEpsTrans(qi, qj))
        {
          dfa_t include_dfa(maximum_refinement(gen, qi_next, qj_next));
          include_dfa.eps_elim();
          ret_dfa.nondet_union(include_dfa);
#if 1
          // Important to keep smaller the automata
          ret_dfa = ret_dfa.makeDFA (alphstart, alphsz).
              minimize (alphstart, alphsz);
              
#endif 
 
        }
      } 
      else 
      {
        // Backwards, an input-consuming transition 
        if(gen.addTrans(qj, gen.word[qj], qi))
        {
          dfa_t include_dfa(maximum_refinement(gen, 
                                               qi_next, 
                                               qj_next));
          include_dfa.eps_elim();
          ret_dfa.nondet_union (include_dfa);
#if 1
          // Important to keep smaller the automata
          ret_dfa = ret_dfa.makeDFA (alphstart, alphsz).
              minimize (alphstart, alphsz);
#endif 
        } 
      }
      gen.pop();
      
      return ret_dfa;
    }
    
  };

} // end namespace covenant

#endif /*__CFG_REFINE_H__*/

