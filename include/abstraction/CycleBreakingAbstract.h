#ifndef __CYCLE_BREAKING_ABSTRACT_H__
#define __CYCLE_BREAKING_ABSTRACT_H__

///////////////////////////////////////////////////////////////////////
// Convert a CFG into a finite automata following the cycle breaking
// heuristics described in the paper "Regular Approximation of CFLs: A
// grammatical view" in "Advances in Probabilistic and Other Parsing
// Technologies volume 16 of Text, Speech, and Language Technology,
// pages 221-241"
///////////////////////////////////////////////////////////////////////
/// A grammar is strongly regular if it does not contain rules of the
/// form A ->* xAy (also called self-embedding rules) where x and y
/// are symbols different from epsilon.  Alternatively, a grammar is
/// strongly regular if for every partition of mutually recursive
/// non-terminals either all productions are left-linear or
/// right-linear.
///
/// The transformation consists of breaking the "unbounded
/// communication" between x and y (called also right and left sides
/// of spines) in order to get rid of self-embedding rules.
///////////////////////////////////////////////////////////////////////

#include <CFG.h>
#include <DFA.h>
#include <abstraction/Abstract.h>
#include <abstraction/StronglyRegularToFA.h>

namespace covenant {

using namespace std;

template<typename EdgeSym>
class CycleBreaking 
{
 public:  
  // PRE: g is normalized.  
  // Rules of the form A->aAc must be converted to A->aBc, B->A (we
  // don't want the left of the rule appearing on the right).
  CFG Approximate(CFG g)
  {
    CFGConnect conn;
    conn.scc(g);
    // cout << "Mutually recursive groups in the grammar: " << conn << endl;
    for (unsigned int i=0; i < conn.nGroups(); i++)
    {
      vector<int>  group = conn.group(i);
      // For efficient membership queries
      set<int>     group_set;
      for(unsigned j=0; j<group.size(); j++)
        group_set.insert(group[j]);
        
      if (!isLinearSCC(g, group, group_set)) 
        cycle_breaking_transf(g, group, group_set);
    }
    return g;
  }
    
 private:

  // A symbol is considered terminal if it is actually terminal or
  // if it is a non-terminal but defined in a different mutually
  // recursive group.
  bool isTerminal (const EdgeSym &s, const set<int> &NonTerminals)
  {
    return (!s.isVar () || NonTerminals.count(s.symID ()) == 0);
  }

  bool isNonTerminal(const EdgeSym &s, const set<int> &NonTerminals)
  {
    return !isTerminal(s,NonTerminals);
  }
    
  // Return true if all symbols from Begin to End are terminal
  bool allTerminal (typename vector<EdgeSym>::const_iterator Begin, 
                    typename vector<EdgeSym>::const_iterator End, 
                    const set<int> &scc)
  {
    for (typename vector<EdgeSym>::const_iterator I = Begin; I != End; ++I)
    {
      if (!isTerminal(*I,scc)) return false;
    }
    return true;
  }

  // A rule is left-linear if it is either A -> A x or A -> x where
  // x is zero or more terminal symbols.
  bool isLLRule (const CFG &g, unsigned int lhs, const vector<EdgeSym> &rhs, 
                 const set<int> &scc)
  {
    if (rhs.empty()) 
      return true;
    if (allTerminal(rhs.begin()+1, rhs.end(), scc))
      return ((rhs[0] == EdgeSym::mkVar(lhs)) || isTerminal(rhs[0],scc));
    else
      return false;
  }

  // A rule is right-linear if it is either A -> x A  or A -> x where
  // x is zero or more terminal symbols.
  bool isRLRule (const CFG &g, unsigned int lhs, const vector<EdgeSym> &rhs, 
                 const set<int> &scc)
  {
    if (rhs.empty()) 
      return true;
      
    if (allTerminal(rhs.begin(), rhs.end()-1, scc))
    {
      return ((rhs[rhs.size()-1] == EdgeSym::mkVar(lhs)) || 
              isTerminal(rhs[rhs.size()-1],scc));
    }
    else
      return false;
  }

  // Return true if all productions in SCC are either left or right
  // linear.
  bool isLinearSCC (const CFG &g, 
                    const vector<int> &group, 
                    const set<int>    &group_set)
  {
    for ( unsigned int vv = 0; vv < group.size(); vv++ )
    {
      for ( unsigned int ri = 0; ri < g.prods[vv].size(); ri++ )
      {
        int r = g.prods[vv][ri].rule;
        vector<EdgeSym> RHS = g.rules[r];
        if (!isLLRule(g, vv, RHS, group_set)) return false;
        if (!isRLRule(g, vv, RHS, group_set)) return false;
      }
    }
    // here all rules are either left-linear or right-linear.
    return true;
  }
    
  // Scan the rule until a non terminal is found.
  // Return true if the whole rule has been scanned.
  bool findNextNonTerminal (const vector<EdgeSym> &r, const set<int> &group, 
                            vector<EdgeSym> &terminal, int & curr)
  {
    unsigned i;
    terminal.clear();
    for(i=curr; i< r.size(); i++)
    {
      if (isTerminal(r[i], group))
        terminal.push_back(r[i]);
      else
        break;
    }
    // mark if the whole rule has been scanned
    if (i < r.size())
    {
      curr = i;
      return false;
    }
    else
      return true;
  }
    
  /*
    It applies the following transformation:

    foreach nonterminal A do
       introduce a new nonterminal A' and add A' -> empty
    foreach prod p: A -> x_0 B1 x_1 B2 ... x_m-1 Bm x_m  (m>=0) do
       replace p with A     -> x_0 B1
       B1'   -> x_1 B2
       ...
       Bm_1' -> x_m-1 Bm
       Bm'   -> x_m A'  
       if A -> a then A -> a A'
  */
  void  cycle_breaking_transf (CFG &g, 
                               const vector<int> &group, 
                               const set<int> &group_set)
  {
    if (group.size() == 1) return;

    boost::unordered_map<int, EdgeSym> new_non_terminals;
    for (unsigned nt=0; nt<group.size(); nt++)
    {
      // add epsilon productions
      EdgeSym A( g.newVar () );
      g.prod(A, Rule::E (g.getTermFactory ()));
      new_non_terminals.insert(make_pair(group[nt], A));
      // LOG ("abstraction" , 
      //      cout << "Added " << A << " -> e " <<  " from " 
      //           << EdgeSym::mkVar(group[nt]) << "\n");
           
    }
    // iterate over each nonterminal in the mutually recursive group 
    for (unsigned nt=0; nt<group.size(); nt++)
    {
      vector<CFG::ProdInfo> pi(g.prods[group[nt]]);
      g.prods[group[nt]].clear();
      // iterate over each rule in the production

      // LOG ("abstraction", 
      //      cout << "Production: "; g.printProduction(cout, group[nt]); cout << endl);

      for( unsigned int ri = 0; ri < pi.size(); ri++ )
      {
        vector<EdgeSym> r(g.rules[pi[ri].rule]);

        // LOG("abstraction", 
        //     cout << "\tRule: "; Rule(r, g.getTermFactory ()); cout << endl);

        cycle_breaking_transf_rule(g, 
                                   group[nt], 
                                   r, 
                                   new_non_terminals, 
                                   group_set);
      }
    }
  }

  // Apply the transfomation
  void cycle_breaking_transf_rule (CFG &g, 
                                   const int lhs, 
                                   const vector<EdgeSym> rule, 
                                   boost::unordered_map<int, EdgeSym> new_non_terminals,
                                   const set<int> group_set)
  {
    int curr_nt=0;
    vector<EdgeSym> terminal_symbols;
    bool finished = findNextNonTerminal (rule, group_set, 
                                         terminal_symbols, curr_nt);
    Rule new_r (g.getTermFactory ());
    for(unsigned i=0; i< terminal_symbols.size(); i++)
      new_r << terminal_symbols[i];
      
    if (finished)
    {
      // the rule has only terminal symbols
      new_r << new_non_terminals[lhs];
      g.prod(EdgeSym::mkVar(lhs), new_r);

      // LOG ("abstraction", 
      //      cout << "\t\trule only with terminal symbols: ";
      //      g.printProduction(cout, lhs); cout << endl);

      return;
    }
    else
    {
      // we found the first nonterminal symbol
      new_r << rule[curr_nt];
      g.prod(EdgeSym::mkVar(lhs), new_r);	

      // LOG ("abstraction" , 
      //      cout << "\t\tAdded rule from first nonterminal: ";
      //      g.printProduction(cout, lhs); cout << endl);

    }
    int prev_nt = curr_nt;
    curr_nt++;
    while (!finished)
    {
      terminal_symbols.clear();
      finished = findNextNonTerminal(rule, group_set, 
                                     terminal_symbols, curr_nt);
      Rule new_rr (g.getTermFactory ());
      for(unsigned i=0; i< terminal_symbols.size(); i++)
        new_rr << terminal_symbols[i];
      if (!finished)
      {
        // intermediate non terminal symbol
        new_rr << rule[curr_nt];
        g.prod(new_non_terminals[rule[prev_nt].symID ()], new_rr);

        // LOG ("abstraction", 
        //      cout << "\t\tintermediate nonterminal symbol: ";
        //      g.printProduction(cout, 
        //                        symID(new_non_terminals[symID(rule[prev_nt])]));
        //      cout << endl);

        prev_nt = curr_nt;
        curr_nt++;
      }
      else
      {
        // last non terminal symbol
        new_rr << new_non_terminals[lhs];
        g.prod(new_non_terminals[rule[prev_nt].symID ()], new_rr);

        // LOG ("abstraction", 
        //      cout << "\t\tlast nonterminal symbol: ";
        //      g.printProduction(cout, 
        //                        symID(new_non_terminals[symID(rule[prev_nt])]));
        //      cout << endl);

        break;
      }
    }
  }
};  // end class CycleBreaking

template< typename EdgeSym>
class CycleBreakingAbstract: public Abstract<EdgeSym> 
{
 public:

  DFA<EdgeSym> do_abstraction(CFG &g, const bool is_regular)
  {
    bool is_abs_done = false;
    DFA<EdgeSym> sfa (g.getTermFactory ());
    if (is_regular)
    {        
      LOG( "verbose", cout << "The CFG is regular thus no need to abstract.\n");

      StronglyRegularToFA<EdgeSym>  conv(g);
      sfa = conv.make_fa();
    }
    else
    {
      LOG( "verbose" , cout << "Applied cycle breaking abstraction.\n");
      is_abs_done = true;
      CycleBreaking<EdgeSym> transformer;
      CFG str_reg_cfg = transformer.Approximate(g);

      LOG ("abstraction" , 
           cout << "Approximated regular grammar:\n" << str_reg_cfg);

      StronglyRegularToFA<EdgeSym> conv(str_reg_cfg);
      sfa = conv.make_fa();
    }

    /// Simplify the automata
    sfa.eps_elim();

    DFA<EdgeSym> dfa     = sfa.makeDFA (g.alphstart, g.alphsz);
    DFA<EdgeSym> min_dfa = dfa.minimize (g.alphstart, g.alphsz);
    bool is_not_empty    = min_dfa.eliminateDeadStates();
    if (!is_not_empty)
    {
      if  (!is_abs_done) 
        throw Exit("One CFG is trivially empty.\n======\nUNSAT\n======");
      else
        throw Exit("One CFG is empty after abstraction.\n======\nUNSAT\n======");
    }
    return min_dfa;
  }
};

} // end namespace covenant

#endif /*__CYCLE_BREAKING_ABSTRACT_H__*/

