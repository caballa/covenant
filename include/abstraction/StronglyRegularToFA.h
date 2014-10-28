#ifndef __STRONGLY_REGULAR_TO_DFA_H__
#define __STRONGLY_REGULAR_TO_DFA_H__

///////////////////////////////////////////////////////////////////////////
/// Convert a strongly regular grammar to a finite automata
///
/// Based on the paper "Practical Experiments with Regular
/// Approximation of Context Free Languages" by M.J. Nederhof.
///////////////////////////////////////////////////////////////////////////

#include <boost/unordered_map.hpp>

#include <Common.h>
#include <DFA.h>
#include <abstraction/Abstract.h>

using namespace std;

namespace covenant  {
    
enum RecursiveType { LEFT, RIGHT, CYCLIC, NON_REC };

template<class EdgeSym>
class StronglyRegularToFA 
{

  typedef boost::unordered_map<unsigned int,unsigned int> SymStateMap;
  typedef boost::unordered_map<unsigned int,unsigned int> SymSCCMap;

  CFG g;                 // g must be a strongly regular grammar
  CFGConnect conn;       // to store the mutually recursive groups of the grammar.
  vector<RecursiveType> recursive;
  SymSCCMap sym_scc_map;
  DFA<EdgeSym> dfa; // the approximation of g.
  bool HasBeenComputed;  // if the approximation was already computed.

 public:  

  StronglyRegularToFA(CFG g_): g(g_), HasBeenComputed(false)
  {
    // The constructor performs three tasks:
    // - compute the strongly connected components of g
    // - compute the map "recursive" which given a SCC returns LEFT,
    //   RIGTH, or CYCLIC, and
    // - precompute some datastructures needed later.
    conn.scc(g);
    // Here to figure out if a mutually recursive group is either:
    // left generating, right generating, or neither.

    vector<bool> left_gen;
    vector<bool> right_gen;
    for (unsigned int i=0; i < conn.nGroups(); i++)
    {
      left_gen.push_back(false);
      right_gen.push_back(false);
    }      

    LOG ("abstraction", 
         cout << "After converted to strongly regular grammar: " 
         << endl << g << endl);

    // FIXME: expensive step just to know which nonterminals generate
    // symbols. Without epsilon-elimination computeGenerating is too
    // conservative.
    CFG g1(g); g1.eliminateEpsilon();
    set<EdgeSym> generatingSymbols = g1.computeGenerating();

    LOG ("abstraction", 
         cout << "Generating symbols: " << generatingSymbols << endl);

    // Traverse all grammar productions
    for ( unsigned int vv = 0; vv < g.prods.size(); vv++ )
    {
      for ( unsigned int ri = 0; ri < g.prods[vv].size(); ri++ )
      {
        int r = g.prods[vv][ri].rule;
        vector<EdgeSym> RHS = g.rules[r];
        // Traverse all recursive groups
        for (unsigned int mi = 0; mi < conn.nGroups() ; mi++){
          vector<EdgeSym> symbols = convert_ints_to_symbols(conn.group(mi));
          // symbols must be sorted
          EdgeSym A = EdgeSym::mkVar(vv);
          if (binary_search(symbols.begin(), symbols.end(), A))
          {
            // mkVar(vv) is a nonterminal from the same SCC
            int pos_first_common_nt = intersect_nth(RHS, symbols);
            if (pos_first_common_nt > -1)
            {
              // EdgeSym B= RHS[pos_first_common_nt];
              // The production rule is of the form: A -> X B Y
              // If  X!=epsilon and Y!=epsilon then marked as left and right
              // If  X!=epsilon but Y==epsilon then marked as left
              // If  X==epsilon but Y!=epsilon then marked as right
              for (int ii=0; (ii<pos_first_common_nt && !left_gen[mi]);ii++)
              {
                if (generatingSymbols.count(RHS[ii]))
                  left_gen[mi]=true; 
              }
              for (int ii=pos_first_common_nt+1;
                   (ii<(int)RHS.size() && !right_gen[mi]);ii++)
              {
                if (generatingSymbols.count(RHS[ii]))
                  right_gen[mi]=true; 
              }
            }
            // This code assumes in the rule A -> X B Y that X and Y
            // are for sure generating but it might not be the case.
            // if (pos_first_common_nt > -1 && 
            //     generatingSymbols.count(RHS[pos_first_common_nt]) > 0){
            //   if (pos_first_common_nt > 0)
            //     left_gen[mi] = true; 
            //   if (pos_first_common_nt < (int) (RHS.size() -1))
            //     right_gen[mi] = true;
            // }                
          }
        }
      }
    }

    LOG ("abstraction", 
         print_left_or_right_generating(left_gen,right_gen));

    // Here we create the "recursive" map. We map each mutually
    // recursive group to either non recursive, left, right, or
    // cyclic.
    for (unsigned int i=0; i < conn.nGroups(); i++)
    {
      if (conn.group(i).size() == 1)
        recursive.push_back(NON_REC);
      else if (!left_gen[i] && right_gen[i])
        recursive.push_back(LEFT);
      else if (left_gen[i] && !right_gen[i])
        recursive.push_back(RIGHT);
      else if (!left_gen[i] && !right_gen[i])
        recursive.push_back(CYCLIC);
      else
      { 
        // Since the grammar is supposed to be strongly regular it
        // cannot be SELF (i.e., both left and right generating).
        throw error("during abstraction from cfg to reg a recursive group is marked as SELF");
      }
    }
      
    LOG ("abstraction", print_recursive());
      
    // reverse map for efficient queries
    // sym_scc_map.reserve(conn.nGroups());

    for (unsigned int i=0; i < conn.nGroups(); i++)
    {
      for (unsigned int j=0; j < conn.group(i).size(); j++)
      {
        SymSCCMap::iterator It = sym_scc_map.find(conn.group(i)[j]);
        if (It == sym_scc_map.end())
          sym_scc_map.insert(make_pair(conn.group(i)[j],i));
        else
          assert((*It).second == i);
      }
    }
  }
    
  // Generate a new state in the DFA and remembers that the state is
  // associated with grammar symbol S.
  inline State MapSymToState(SymStateMap &M, EdgeSym S)
  {
    SymStateMap::iterator It = M.find(S.symID ());
    if (It != M.end())
      return mkState((*It).second);
    else
    {
      State q = dfa.state(dfa.nStates());
      M.insert(make_pair(S.symID (), q.id));
      return q;
    }
  }
    
  // Convert the CFG into a strongly regular language in a recursive
  // manner. It is a direct implementation of the algorithm
  // described in the paper "Practical Experiments with Regular
  // Approximation of Context Free Languages" by M.J. Nederhof.
  //
  // FIXME: this version can have an exponential blowup. In the
  // abovementioned paper there is a solution to avoid the blowup.
  void make_fa(const State q_0, 
               const vector<EdgeSym> alpha, 
               const State q_1, 
               SymStateMap &SymToState) 
  {

    LOG ("abstraction", 
         cout << "reading word " << alpha 
         << " with size " << alpha.size() << endl);

    if (alpha.size() == 0)
    {

      LOG ("abstraction", 
           cout << "[make_fa] reading epsilon\n");

      dfa.eps_transition(q_0,q_1);
      return;
    }
    else if (alpha.size () == 1 && alpha[0].isTerm ())
    {

      LOG ("abstraction", 
           cout << "[make_fa] (case A->a) reading " << alpha << endl);

      dfa.transition(q_0, alpha[0], q_1);
      return;
    }
    else if (alpha.size() >= 2)
    {
      vector<EdgeSym> beta(alpha);
      remove(beta,0);
      State q_2 = dfa.state(dfa.nStates());
      vector<EdgeSym> gamma;
      gamma.push_back(alpha[0]);

      LOG ("abstraction",
           cout << "[make_fa] (case A-> BC) reading " << gamma << " " 
           << beta << endl);

      make_fa(q_0, gamma, q_2, SymToState);
      make_fa(q_2, beta , q_1, SymToState);
      return;
    }
    else
    {
      //assert(alpha.size() == 1 && isVar(alpha[0]));

      LOG ("abstraction", 
           cout << "[make_fa] (case A->B) reading " << alpha << endl);

      EdgeSym A =  alpha[0];	
      SymSCCMap::iterator It = sym_scc_map.find(A.symID ());
      if (recursive[(*It).second] == NON_REC)
      {

        LOG ("abstraction", 
             cout << "[make_fa] A is not recursive\n ");

        for ( unsigned int ri = 0; ri < g.prods[A.symID()].size(); ri++ )
        {
          int r = g.prods[A.symID ()][ri].rule;
          vector<EdgeSym> RHS = g.rules[r];
          make_fa(q_0, RHS, q_1, SymToState);
        }
        return;
      }
        
      State q_A = MapSymToState(SymToState, A);
      if (recursive[(*It).second] == LEFT)
      {

        LOG ("abstraction", 
             cout << "[make_fa] A is marked as LEFT\n");

        // GroupSyms must be sorted
        vector<EdgeSym> GroupSyms = convert_ints_to_symbols(conn.group((*It).second));
        for (unsigned int gi =0; gi < GroupSyms.size() ; gi ++)
        {
          EdgeSym C = GroupSyms[gi];
          State q_C = MapSymToState(SymToState, C);
          for ( unsigned int ri = 0; ri < g.prods[C.symID ()].size(); ri++ )
          {
            int r = g.prods[C.symID ()][ri].rule;
            vector<EdgeSym> RHS = g.rules[r];
              
            if (!intersect_unorder_ord(RHS, GroupSyms))
              make_fa(q_0, RHS, q_C, SymToState);
            else{
              if (binary_search(GroupSyms.begin(),GroupSyms.end(), RHS[0]) &&
                  (!intersect_unorder_ord(++RHS.begin(), RHS.end(), GroupSyms)))
              {
                EdgeSym D = RHS[0];
                vector<EdgeSym> beta(RHS);
                remove(beta,0);
                State q_D = MapSymToState(SymToState, D);
                make_fa(q_D, beta, q_C, SymToState);
              }		    
            }
          }
        }
        dfa.eps_transition(q_A,q_1);
      } // end LEFT
      else
      {
        /////
        // RIGHT or CYCLIC: code almost identical
        /////

        LOG ("abstraction", 
             cout << "[make_fa] A is marked as RIGHT/CYCLIC\n"); 

        vector<EdgeSym> GroupSyms = convert_ints_to_symbols(conn.group((*It).second));
        for (unsigned int gi =0; gi < GroupSyms.size() ; gi ++)
        {
          EdgeSym C = GroupSyms[gi];
          State q_C = MapSymToState(SymToState, C);
          for ( unsigned int ri = 0; ri < g.prods[C.symID ()].size(); ri++ )
          {
            int r = g.prods[C.symID ()][ri].rule;
            vector<EdgeSym> RHS = g.rules[r];

            LOG ("abstraction", cout << RHS << "   " );

            if (!intersect_unorder_ord(RHS, GroupSyms))
              make_fa(q_C, RHS, q_1, SymToState);
            else
            {
              if (binary_search(GroupSyms.begin(),GroupSyms.end(), RHS[RHS.size()-1]) &&
                  (!intersect_unorder_ord(RHS.begin(), RHS.end()-1, GroupSyms)))
              {
                EdgeSym D = RHS[RHS.size()-1];
                vector<EdgeSym> beta(RHS);
                remove(beta,RHS.size()-1);
                State q_D = MapSymToState(SymToState, D);
                make_fa(q_C, beta, q_D,  SymToState);
              }		    
            }
          }
        }
        dfa.eps_transition(q_0,q_A);
      } // end RIGHT or CYCLIC	
      return;
    }
  }

  DFA<EdgeSym> make_fa()
  {
    if (HasBeenComputed) 
      return dfa;
    else{
      State start, final;
      start = dfa.state(0);
      final = dfa.state(1);           
      dfa.setStart(start);
      dfa.accept(final);
      vector<EdgeSym> alpha;
      alpha.push_back(g.startSym());
      SymStateMap SymToState;
      // SymToState.reserve(g.nVars);
      make_fa(start, alpha, final, SymToState);
      HasBeenComputed=true;
      return dfa;
    }
  }

 private:

  // Convert a vector of integers into a vector of symbols This
  // transformation does not affect the order.
  vector<EdgeSym> convert_ints_to_symbols(vector<int> symbols_ids)
  {
    vector<EdgeSym> symbols;
    transform(symbols_ids.begin(), symbols_ids.end(), 
              back_inserter(symbols), EdgeSym::mkVar);
    return symbols;
  }
    
  // Return -1 if RHS and sorted_seq has no element in
  // common. Otherwise, it returns the position of the first common
  // element.  Complexity is O(mlog(n)) here m=|RHS| and n=|sorted_seq|.
  int intersect_nth(const vector<EdgeSym> &RHS, const vector<EdgeSym> &sorted_seq)
  {
    for (unsigned i=0; i < RHS.size() ; i++)
    {
      if (binary_search(sorted_seq.begin(), sorted_seq.end(), RHS[i]))
      {
        return i;
      }
    }
    return -1;
  }

  // For logging
  void print_recursive()
  {
    for (unsigned int i=0; i < conn.nGroups(); i++)
    {
      cout << "The mutually recursive group { ";
      for (unsigned int j=0 ; j < conn.group(i).size() ; j++)
      {
        cout << EdgeSym::mkVar(conn.group(i)[j]) << " ";
      }
      cout << "} is: ";
      switch(recursive[i])
      {
        case NON_REC:
          cout << "NOT RECURSIVE" << endl;
          break;
        case LEFT:
          cout << "LEFT" << endl;
          break;
        case RIGHT:
          cout << "RIGHT" << endl;
          break;
        case CYCLIC:
          cout << "CYCLIC" << endl;
          break;
        default: ;;	 
      }
    }
  }
    
  // For logging
  void print_left_or_right_generating(const vector<bool> &left_gen,
                                      const vector<bool> &right_gen)
  {
    for (unsigned int i=0; i < conn.nGroups(); i++)
    {
      cout << "The mutually recursive group { ";
      for (unsigned int j=0 ; j < conn.group(i).size() ; j++)
      {
        cout << EdgeSym::mkVar(conn.group(i)[j]) << " ";
      }
      cout << "} is marked as: ";
        
      if (left_gen[i]) 
        cout << "left generating " << endl;
      if (right_gen[i]) 
        cout << "right generating " << endl;
      if (!left_gen[i] &&  !right_gen[i])
        cout << "only terminals in the rule " << endl;
    }
  }

};  

} // end namespace covenant

#endif /* __STRONGLY_REGULAR_TO_DFA_H__ */
