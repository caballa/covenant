#ifndef __CFG_COND_EPS_GEN_H__
#define __CFG_COND_EPS_GEN_H__
/////
// Code for Epsilon-generalization with respect to a CFG L.
// Performs roughly the same task as pre*; just
// a slimmed-down implementation.
/////
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>

#include <Common.h>
#include <CFG.h>
#include <DFA.h>

#define EPSSYM (-1)

namespace covenant {

  using namespace std;

  // Precomputed data-structures for
  // doing efficient intersection
  // tests on CFGs.

  // Possibly worth trying to compact the CFG.
  template<typename EdgeSym>
  class CFGDigest {
  public:
    typedef map_set<int, int, boost::hash<int>, std::equal_to<int> > sing_map;
    typedef pair<int, int> spair;
    // Precondition -- the CFG is in near-CNF.
    // So, only singleton productions may use terminals,
    // and productions have degree at most 2.
    // It differs from CNF in that
    //  - empty productions are permitted
    //  - singleton non-terminal productions are permitted
    CFGDigest(CFG& g)
      : nprods(g.prods.size()),
        startSym(g.startSym().symID ()),
        extend_sing(nprods, vector<int>()),
        extend_fwd(nprods, vector<spair>()),
        extend_rev(nprods, vector<spair>())
    {
      for(int pi = 0; pi < nprods; pi++)
      {
        vector<CFG::ProdInfo>& ps(g.prods[pi]);
        for(unsigned int ri = 0; ri < ps.size(); ri++)
        {
          vector<EdgeSym>& r(g.rules[ps[ri].rule]);

          if (r.size() == 0)
          {
            // Nullable production.
            nullable_prods.push_back(pi);
          } 
          else if (r.size() == 1) {
            // Singleton production.
            EdgeSym s(r[0]);
            if (s.isTerm ())
            {
              term_prods[s.symID ()].insert(pi);
            } 
            else {
              extend_sing[s.symID ()].push_back(pi);
            }
          } 
          else {
            assert(r.size() == 2);
            int b = r[0].symID ();
            int c = r[1].symID ();
            extend_fwd[b].push_back(spair(c, pi));
            extend_rev[c].push_back(spair(b, pi));
          }
        }
      }
    }
   
    int nprods;
    int startSym;
    vector<int> nullable_prods;
    sing_map term_prods;

    vector< vector<int> > extend_sing;
    vector< vector<spair> > extend_fwd;
    vector< vector<spair> > extend_rev;
  }; // end class CFGDigest

  template<typename EdgeSym>
  class CondEpsGen {
    
    typedef struct 
    {
      int start;
      int end;
      int s;
    } prod_info;

    prod_info mk_prod(int start, int s, int end) 
    {
      prod_info p = { start, end, s };
      return p;
    }

  public:

    typedef DFA<EdgeSym> dfa_t;
    typedef CFGDigest<EdgeSym> cfg_digest_t;

    CondEpsGen(cfg_digest_t& _digest, const vector<int>& _word)
        : digest(_digest), word(_word),
          nstates(1+word.size()), nprods(digest.nprods),
          produce_fwd(vector< vector<int> >(nstates*nprods, vector<int>())),
          produce_rev(vector< vector<int> >(nstates*nprods, vector<int>())),
          produced(vector<bool>(nstates*nstates*nprods, false)),
          //      produced((bool*) calloc(nstates*nstates*nprods, sizeof(bool))),
          head(0),
          eps_preds(nstates, vector<int>()),
          eps_succs(nstates, vector<int>())
    {
      // Assumption -- g is semi-normalized.
      // Single productions may be nullable, but
      // terminals may only occur in singleton productions; this means
      // we only need |P||Q|^2 space, rather than (|P|+|Sigma|)|Q|^2.
      
      // Initialize the leaf transitions.  
      for(unsigned int ci = 0; ci < word.size(); ci++)
      {
        set<int>& sing(digest.term_prods[word[ci]]);
        
        set<int>::iterator it;
        for(it = sing.begin(); it != sing.end(); ++it)
        {
          addProd(ci, *it, ci+1);
        }
      }
      
      // Add empty productions
      for(unsigned int ci = 0; ci <= word.size(); ci++)
      {
        for(unsigned int pi = 0; pi < digest.nullable_prods.size(); pi++)
        {
          addProd(ci, digest.nullable_prods[pi], ci);
        }
      }
      // Since we know we only recognize w, we can either
      // use the standard propagation approach, or use
      // a standard(-ish) CYK implementation.
      
      // For now, we'll use the worklist approach, since
      // it's already implemented.
      update_prods(); 
    }

    ~CondEpsGen() { }

    // Returns true if the language is empty
    bool addEpsTrans(int qi, int qj)
    {
      assert(qi < qj);
      
      // FIXME: Check that the epsilon transition doesn't
      // already (transitively) exist.

      // Record the decision
      decisions.push_back(mk_prod(qi, EPSSYM, qj));

      // Add it to the set of transitions
      eps_preds[qj].push_back(qi);
      eps_succs[qi].push_back(qj);
      
      // Extend any {qh, A, qi} to {qh, A, qj}
      for(int p = 0; p < nprods; p++)
      {
        vector<int>& rev(prod_rev(qi, p));
        for(unsigned int pi = 0; pi < rev.size(); pi++)
        {
          int qh = rev[pi];
          addProd(qh, p, qj);
        }
      }
      
      // Now extend any {qj, A, qk} to {qi, A, qk}
      for(int p = 0; p < nprods; p++)
      {
        vector<int>& fwd(prod_fwd(qj, p));
        for(unsigned int pi = 0; pi < fwd.size(); pi++)
        {
          int qk = fwd[pi];
          addProd(qi, p, qk); 
        }
      }
      
      return update_prods();
    }
    
    bool addTrans(int qi, int c, int qj)
    {
      // Update the automaton with a new transition.
      // Not presently checking if p[qi -> qj] alreacy
      // exists.
      decisions.push_back(mk_prod(qi, c, qj));
      
      set<int>& sing(digest.term_prods[c]);
      
      set<int>::iterator it;
      for(it = sing.begin(); it != sing.end(); ++it)
      {
        addProd(qi, *it, qj);
      }
      return update_prods();
    }
    
    // Convert an augmented-word into an automaton.
    dfa_t emit_dfa(void){
      dfa_t dfa;
      
      for(unsigned i=0; i < word.size(); i++){
        // This is gonna at least double the size of the alphabet. O__O
        dfa.transition(dfa.state(i), EdgeSym::mkTerm(word[i]), dfa.state(i+1));
      }
      dfa.setStart(dfa.state(0));
      dfa.accept(dfa.state(word.size()));
      
      // Add augmentations
      for(unsigned int di = 0; di < decisions.size(); di++)
      {
        prod_info& d(decisions[di]); 
        if(d.s == EPSSYM)
        {
          dfa.eps_transition(dfa.state(d.start), dfa.state(d.end));
        } else {
          dfa.transition(dfa.state(d.start), EdgeSym::mkTerm(d.s), dfa.state(d.end));
        }
      }
      
      return dfa;
    }

    void push(void)
    {
      lim.push_back(added.size());
      declim.push_back(decisions.size());
    }

    void pop(void)
    {
      assert(lim.size() > 0);
      int limit = lim.back();
      assert(limit >= 0);
      lim.pop_back(); 

      head = added.size();
      while(limit < head)
      {
        prod_info p = added[--head];

        prod(p.start, p.s, p.end) = false;
        prod_fwd(p.start, p.s).pop_back();
        prod_rev(p.end, p.s).pop_back(); 
      }
      added.resize(limit);
      assert(head >= 0);

      int dlim = declim.back();
      declim.pop_back();
      decisions.resize(dlim);
    }

    bool is_empty(void)
    {
      // Check if S is generated by (q0 -> qN)
      return !prod(0, digest.startSym, nstates-1);
    }

    int nStates(void) { return nstates; }

   protected:

    void addProd(int qi, int p, int qj)
    {
      if(!prod(qi, p, qj))
      {
        prod(qi, p, qj) = true;
        prod_fwd(qi, p).push_back(qj);
        prod_rev(qj, p).push_back(qi);
        added.push_back(mk_prod(qi, p, qj));
      }
    }
    
    // Returns true if the language is empty
    bool update_prods(void)
    {
      prod_info p;
      assert(head >= 0);
      while(head < (int) added.size())
      {
        // Pick the next unprocessed production
        p = added[head++]; 
        
        int qi = p.start;
        int qj = p.end;
        int s = p.s;
        
        // Check for new productions
        
        // Epsilon-transitions
        for(unsigned int pred = 0; pred < eps_preds[qi].size(); pred++)
        {
          addProd(eps_preds[qi][pred], s, qj);
        }
        for(unsigned int succ = 0; succ < eps_succs[qj].size(); succ++)
        {
          addProd(qi, s, eps_succs[qj][succ]);
        }
        
        // Singleton extensions
        for(unsigned int pi = 0; pi < digest.extend_sing[s].size(); pi++)
          addProd(qi, digest.extend_sing[s][pi], qj);
        
        // Forward: B[qi -> qj] C[qj -> qk] |- A [qi -> qk]
        for(unsigned int fi = 0; fi < digest.extend_fwd[s].size(); fi++)
        {
          // Satisfying (C, A) pair.
          typename cfg_digest_t::spair ca = digest.extend_fwd[s][fi];
          
          vector<int>& ks(prod_fwd(qj, ca.first)); 
          for(unsigned int kk = 0; kk < ks.size(); kk++)
            addProd(qi, ca.second, ks[kk]);
        }
        // Backward
        for(unsigned int bi = 0; bi < digest.extend_rev[s].size(); bi++)
        {
          typename cfg_digest_t::spair ba = digest.extend_rev[s][bi];
          
          vector<int> hs(prod_rev(qi, ba.first));
          for(unsigned int hh = 0; hh < hs.size(); hh++)
            addProd(hs[hh], ba.second, qj);
        }
      }
      
      return is_empty();
    }

    inline vector<bool>::reference prod(int qi, int p, int qj)
    {
      return produced[nstates*(nstates*p + qi) + qj];
    }
    /*
    inline bool& prod(int qi, int p, int qj)
    {
      assert(qi < nstates);
      assert(qj < nstates);
      assert(p < nprods);
      //  return produced[nstates*(nstates*p + qi) + qj];
      return produced[nstates*nstates*p + nstates*qi + qj];
    }
    */

    inline vector<int>& prod_fwd(int qi, int p)
    {
      return produce_fwd[nstates*p + qi];
    }

    inline vector<int>& prod_rev(int qj, int p)
    {
      return produce_rev[nstates*p + qj];
    }

    cfg_digest_t& digest;
 public:
    const vector<int> word;
 protected:
    
    int nstates;
    int nprods;

    // Precomputed extensions
    // extend_fwd: A +-> (B, C) iff C -> A B in g.

    // Lookup tables
    // produce_fwd: q_i, A +-> q_j iff g |- A over (q_i, q_j).
    // produce_rev: q_j, A +-> q_i iff g |- A over (q_i, q_j).
    // produced: (q_i, A, q_j) = T iff g |- A over (q_i, q_j).
    vector< vector<int> > produce_fwd;
    vector< vector<int> > produce_rev;
    vector<bool> produced;
    // bool* produced;

    // Trailing
    int head; // Next production to process
    vector<int> lim;
    vector<prod_info> added; // Productions that were added.

    vector<int> declim;
    vector<prod_info> decisions;

    vector< vector<int> > eps_preds;
    vector< vector<int> > eps_succs;
  }; // end CondEpsGen

}
#endif /* __CFG_COND_EPS_GEN_H__ */
