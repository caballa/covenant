#ifndef __DFA_H__
#define __DFA_H__

#include <cassert>
#include <iostream>
#include <vector>
#include <set>

#include <boost/optional.hpp>

#include <adt/SparseSet.h>
#include <adt/UnionFind.h>

//=======================================================================
// Finite automaton implementation
//=======================================================================
// While the code says DFA, it's essentially an NFA implementation,
// since nothing is preventing there from being overlapping transitions
// from a state.
//=======================================================================

namespace covenant {

  namespace DFA_impl{

     // Data structure for hashing SFA states.
     template<typename Container>
     class StateTable 
     {
       typedef struct 
       {
         unsigned int sz;
         unsigned int data[1];
       } StateEntry;
       
       struct HashEntry {
         std::size_t operator()(const StateEntry* x) const{
           std::size_t res = 0;
           for(unsigned int i = 0; i < x->sz; i++){
             boost::hash<const unsigned int> hasher;
             boost::hash_combine(res, hasher(x->data[i]));
           }
           return res;
         }
       };
       
       struct EqEntry {
         bool operator()(const StateEntry* x, const StateEntry* y) const {
           if(x->sz != y->sz)
             return false;
           for(unsigned int ii = 0; ii < x->sz; ii++)
             if(x->data[ii] != y->data[ii]) return false;
           return true;
         }
       };
       
       typedef boost::unordered_map<const StateEntry*, int, HashEntry, EqEntry> StateCache;
       
       int tempsz;
       StateEntry* temp;
       vector<StateEntry*> entries;
       StateCache cache; 
       
      public:
       
       StateTable(void): tempsz(1)
       {
         temp = ((StateEntry*) malloc(sizeof(StateEntry) + (2*tempsz)*sizeof(unsigned int)));
       }
       
       ~StateTable()
       {
         // Free allocated memory.
         free(temp);
         for(unsigned int ii = 0; ii < entries.size(); ii++)
           free(entries[ii]);
       }
       
       // Return the value of the key if it is already in the
       // table. Otherwise, it returns -1.
       int find(Container & es, bool accept)
       {
         // Ensure there is enough space in temp for the edges
         if(tempsz < (int) es.size())
         {
           while(tempsz < (int) es.size())
             tempsz *= 2;
           temp = ((StateEntry*) realloc(temp, sizeof(StateEntry) + 
                                         (2*tempsz)*sizeof(unsigned int)));
         }
         // Load the edges into temp.
         temp->sz = es.size() + 1;
         temp->data[0] = accept;
         unsigned int ii=0;
         for(typename Container::iterator it = es.begin(), et = es.end(); it!=et; ++it)
         {	
           temp->data[ii+1] = *it;
           ii++;
         }
         
         // Look for the temporary element in the map.
         typename StateCache::iterator it(cache.find(temp));
         if (it != cache.end())
           return (*it).second;
         else
           return -1;
       }
       
       // Check if the given sequence of type Container of states is
       // already in the table.  If so, return the corresponding id;
       // otherwise, insert and return base.
       unsigned int insert(Container & es, bool accept, unsigned int base){
         // Ensure there is enough space in temp for the edges
         if(tempsz < (int) es.size())
         {
           while(tempsz < (int) es.size())
           {
             tempsz *= 2;
           }
           temp = ((StateEntry*) realloc(temp, sizeof(StateEntry) + 
                                         (2*tempsz)*sizeof(unsigned int)));
         }
         // Load the edges into temp.
         temp->sz = es.size() + 1;
         temp->data[0] = accept;
         unsigned int ii=0;
         for(typename Container::iterator it = es.begin(), et = es.end(); it!=et; ++it)
         {
           temp->data[ii+1] = *it;
           ii++;
         }
         
         // Look for the temporary element in the map.
         typename StateCache::iterator it(cache.find(temp));
         if(it != cache.end())
           return (*it).second;
         // Look for the temporary element in the map.
         if(it != cache.end())
           return (*it).second;
         else 
         {
           StateEntry* entry = ((StateEntry*) malloc(sizeof(StateEntry) + 
                                                     (2*es.size())*sizeof(unsigned int))); 
           entries.push_back(entry);
           entry->sz = es.size() + 1; 
           entry->data[0] = accept;
           unsigned int ii=0;
           for(typename Container::iterator it = es.begin(), et = es.end(); it!=et; ++it)
           {
             entry->data[ii+1] = *it;
             ii++;
           }
           cache[entry] = base;
           return base;
         }
       }
       
       void clear()
       {
         cache.clear();
         for(unsigned int ii = 0; ii < entries.size(); ii++)
         {
           delete entries[ii];
         }
         entries.clear();
       }     
     }; // end class StateTable
  
     template<typename V>
     class distinctTable 
     {
       vector<vector<V> > M;

      public:
       distinctTable(unsigned n, V init_val): M(n) 
       {
         for (unsigned i = 0; i < n; ++i)
           M[i] = vector<V>(n);
         
         for (unsigned i = 0; i < n; ++i)
         {
           for (unsigned j = 0; j < n; ++j)
             M[i][j]=init_val;
         }
       }
       vector<V>& operator[](int i) { return M[i]; }
       unsigned size() const { return M.size(); }
     };
  
     inline ostream& operator<<(ostream& o, distinctTable<bool> T) 
     {
       for(unsigned int i=0; i < T.size(); i++)
       {
         for(unsigned int j=T.size()-1; j > i; j--)
         {
           o << "[" << i << "," << j << "]";
           if (T[i][j]) 
             o << "!= "; 
           else 
             o << "= ";
         }
         o << endl;
       }  
       return o;
     }

  } // end namespace DFA_impl


  // A light wrapper for states, to help avoid accidentally using a state
  // that hasn't been allocated.
  struct State 
  {
    unsigned int id;

    friend State mkState(int id);
    
    bool operator == (State p) const { return id == p.id; }
    bool operator != (State p) const { return id != p.id; }
    bool operator < (State p)  const  { return id < p.id; }
  };

  inline State mkState(int id) { State p; p.id = id; return p; }

  inline std::ostream& operator<<(std::ostream& o, State s) 
  {
    o << s.id;
    return o;
  }

  // A symbolic automaton, parameterized by transition-relation
  // representation.
  template<class V>
  class DFA 
  {

    typedef boost::unordered_map<int , set<unsigned int> > ReverseStateTable;
    typedef boost::unordered_map<unsigned int, unsigned int> MapStates;

  public:

    typedef DFA<V> dfa_t;

    typedef struct {
      V val;
      unsigned int dest;
    } Edge;
    
    Edge mk_edge(const V& val, unsigned int dest) 
    {
      Edge e = { val, dest };
      return e;
    }

    int start;
    std::vector< std::vector<Edge> > trans; // The state transitions.
    std::vector< std::vector<int> > eps;    // Epsilon transitions
    std::vector<bool> accepts;

    DFA() : start( 0 ) { }
    
    void clear()
    {
      trans.clear();
      eps.clear();
      accepts.clear();
      start = 0;
    }

    State state(int i) 
    {
      while((int) trans.size() <= i)
      {
        trans.push_back(std::vector<Edge>());
        eps.push_back(std::vector<int>());
        accepts.push_back(false);
      }
      return mkState(i);
    }

    void setStart(State s)
    {
      start = s.id;
    }

    // Mark a state as accepting.
    // Note that accepts is an array of bool flags, not
    // a list of accepting states.
    void accept(State s)
    {
      assert(s.id < trans.size());
      accepts[s.id] = true;
    }
    
    State startState(void) const 
    {
      return mkState(start);
    }
    
    // (s, val) -> r
    void transition(State s, const V& val, State r) 
    {
      while( s.id >= trans.size() )
        trans.push_back( std::vector<Edge>() );
      
      trans[s.id].push_back( mk_edge(val, r.id) );
    }

    // Add an epsilon transition
    void eps_transition(State s, State r)
    {
      eps[s.id].push_back(r.id);
    }

    // Return true if the automata has already the transition 
    inline bool has_transition(State src, State dest, V val )
    {
      assert( src.id < (unsigned) nStates());
      std::vector<Edge> &ts( trans[src.id] );
      for(unsigned int j=0; j < ts.size(); j++){
        if (ts[j].dest == dest.id && ts[j].val == val)
          return true;
      }
      return false;
    }
    
    std::ostream& write(std::ostream &o)
    {
      o << "[" << start << " => ";
      bool first = true;
      for(unsigned int ai = 0; ai < accepts.size(); ai++)
      {
        if(accepts[ai])
        {
          if(first)
            first = false;
          else
            o << ", ";
          o << ai;
        }
      }
      o << "]" << std::endl;
      for( unsigned int ss = 0; ss < trans.size(); ss++ )
      {
        std::vector<Edge>& ts( trans[ss] );
        o << ss << " [";

        for( unsigned int ti = 0; ti < ts.size(); ti++ )
        {
          if( ti != 0 )
            o << ", ";
	  o << ts[ti].val << " -> " << ts[ti].dest;
        }
        o << "]";
        
        std::vector<int>& es( eps[ss] );
        if(es.size() > 0)
        {
	  o << "[e -> ";
          o << es[0];
          for(unsigned int ei = 1; ei < es.size(); ei++)
            o << ", " << es[ei];
          o << "]";
        }
        o << std::endl;
      }
      return o;
    }

    // Print the automata in dot format
    std::ostream& print_dot(std::ostream& o, std::string label) 
    {

      static unsigned int sfa_key=0;
      sfa_key++;
      o << "digraph " << "sfa_" << sfa_key << "{" << std::endl;
      o << "label = \"" << label << "\";" << std::endl;
      o << "rankdir=LR;" << std::endl;
      o << "node [shape = point ] qi;" << std::endl;
      for (int i=0; i < nStates() ; i++)
      {
        if (accepts[i])
          o << "node [shape = doublecircle ] " << "q" << i << ";" << std::endl;
        else
          o << "node [shape = circle ] " << "q" << i << ";" << std::endl;
      }
      o << "qi" << " -> " << "q" << startState().id << std::endl;
      for( unsigned int ss = 0; ss < trans.size(); ss++ )
      {
        for( unsigned int ti = 0; ti < trans[ss].size(); ti++ )
        {
          o << "q" << ss << " -> " << "q" << trans[ss][ti].dest 
            << " [label=\"" << trans[ss][ti].val << "\"];" << std::endl;
        }
        for(unsigned int ei = 0; ei < eps[ss].size(); ei++)
        {
          o << "q" << ss << " -> " << "q" << eps[ss][ei] 
            << " [label=e];" << std::endl;
        }
      }
      o << "}" << std::endl;
      return o;
    }

    // Compute the epsilon closure 
    void eps_close(int state)
    {
      SparseSet<> reach(nStates());
      reach.insert(state);
      for(unsigned int ei = 0; ei < eps[state].size(); ei++)
      {
        // Added the conditional test after observing some problems
        // running the CFG solver.
        if(!reach.elem(eps[state][ei]))
          reach.insert(eps[state][ei]);
      }

      for(unsigned int ii = 1; ii < reach.size(); ii++)
      {
        int es = reach[ii];
        for(unsigned int ei = 0; ei < eps[es].size(); ei++)
        {
          // Not yet processed
          if(!reach.elem(eps[es][ei]))
          { 
            reach.insert(eps[es][ei]);
            eps[state].push_back(eps[es][ei]);
          }
        }
      }
    }

    // Eliminates all epsilon transitions from the automaton.
    // At this stage, doesn't attempt to eliminate redundant or overlapping
    // transitions. 
    void eps_elim(void)
    {
      // Compute the transitive closure of epsilon transitions
      for(int ii = 0; ii < nStates(); ii++)
      {
        eps_close(ii);
      }

      // Determine the set of edges introduced by the epsilon transitions
      // Why do it this way? So we don't include multiple copies when there is a chain
      // of epsilon transitions.
      std::vector< std::vector<Edge> > eps_induced;
      for(int ii = 0; ii < nStates(); ii++)
      {
        eps_induced.push_back(std::vector<Edge>());
        for(unsigned int ei = 0; ei < eps[ii].size(); ei++)
        {
          // Update accept states
          if(accepts[eps[ii][ei]])
            accepts[ii] = true;

          // Add all the epsilon-reachable transitions
          eps_induced[ii].insert(eps_induced[ii].end(),
                                 trans[eps[ii][ei]].begin(), 
                                 trans[eps[ii][ei]].end());
        }
      }

      // Remove the epsilon transitions, and add the induced transitions
      for(int ii = 0; ii < nStates(); ii++)
      {
        eps[ii].clear();
        trans[ii].insert(trans[ii].end(), 
                         eps_induced[ii].begin(), 
                         eps_induced[ii].end());
      }
    }

    int nStates(void) const { return trans.size(); }

    // if changes[i] is false then the i-th state and all the incoming
    // and outgoing edges to/from i-th state are removed.
    // Return true if there is still some accepting state.
    inline dfa_t resize(const vector<bool>  &changes, bool &is_accepting)
    {
      MapStates mapping; // mapping to relabel states old -> new
      unsigned int nstates=0;
      int newstart=-1;

      if (this->start < 0) 
        throw error("DFA does not have an initial state");

      for(unsigned int i=0; i< (unsigned int) this->nStates(); i++)
      {
        if (changes[i]){
          if (i == (unsigned int) this->start) 
            newstart = nstates;
          mapping.insert(make_pair(i, nstates));
          nstates++;
        }
      }

      if (newstart < 0) 
        throw error("DFA does not have an initial state");

      dfa_t resized_fa;    

      // Create new states
      for(unsigned int i=0; i<nstates; i++)
      {
        resized_fa.state(i);
      }

      // Marking initial state
      resized_fa.setStart(mkState((unsigned int) newstart));
      // Copy delta transitions
      for(unsigned int i=0; i< (unsigned int) this->nStates(); i++)
      {
        for(unsigned int j=0; j< this->trans[i].size(); j++)
        {
          if (changes[i] && changes[this->trans[i][j].dest])
          {
            resized_fa.transition(mkState(mapping[i]), 
                                  this->trans[i][j].val,
                                  mkState(mapping[this->trans[i][j].dest]));
          }
        }
      }
      // Marking final states
      is_accepting=false;
      for(unsigned int i=0; i< (unsigned int) this->nStates(); i++){
        if (changes[i] && this->accepts[i]){
          resized_fa.accept(mkState(mapping[i]));
          is_accepting=true;
        }
      }
      return resized_fa;
    }
    
    // Pre: *this has no epsilon transitions
    // Return true iff the language of *this is empty.
    inline bool empty() const
    {
      return !is_generating_state(mkState(this->start));
    }

    // Pre : *this does not have epsilon transitions.
    // Post: eliminate dead states from *this and return true iff is
    //       there is some accepting path.
    bool  eliminateDeadStates()
    {
      // initialize vector
      vector<bool> live_states;
      live_states.reserve(this->nStates());
      for(int i=0; i< this->nStates(); i++)
      {
        live_states.push_back(false);
      }
      bool change=true;
      while (change)
      {
        change=false;
        set<State> visited;
        vector<bool> old_live_states(live_states);
        mark_generating_states(this->startState(), visited, live_states);
        for(unsigned i=0; i<live_states.size();i++)
        {
          if (live_states[i] != old_live_states[i])
            change=true;
        }
      }
      
      if (!live_states[this->startState().id])
      {
        this->clear();
        return false;
      }
      else
      {
        bool is_accepting;
        *this = resize(live_states, is_accepting);
        return is_accepting;
      }
    }

#ifdef SANITY_CHECKS
    void CheckNoEpsilonTrans() const
    {
      for(int ii = 0; ii < this->nStates(); ii++)
      {
        if (!this->eps[ii].empty())
          throw error("DFA cannot have epsilon transitions");
      }
    }
#endif 

    // Convert a NFDA to FDA.
    // Pre: *this does not have epsilon transitions
    //       alphstart and alphsz determine Sigma 
    dfa_t makeDFA(unsigned alphastart, unsigned alphasz) const
    {

#ifdef SANITY_CHECKS    
      CheckNoEpsilonTrans();
#endif 
      
      LOG ("dfa-determinisation", 
           cout << "Converting from NDFA to DFA ... \n");
      
      // The deterministic dfa
      dfa_t dfa;

      // Initialise bookkeeping.
      unsigned int ID = 0;
      DFA_impl::StateTable<set<unsigned int> > table;
      ReverseStateTable rev_table;
      table.clear();
      rev_table.clear();

      // Special "sink" state
      State sink = dfa.state(ID++);
      set<unsigned int> ss;
      ss.insert(this->startState().id);    
      table.insert(ss, this->accepts[this->startState().id], ID);
      rev_table.insert(make_pair(ID,ss));
      vector<unsigned int> worklist;
      worklist.push_back(ID);
      State start = dfa.state(ID++);
      dfa.setStart(start);
      
      while (!worklist.empty())
      {
        unsigned int ss_ID = worklist.back();
        set<unsigned int> ss = rev_table[ss_ID];
        worklist.pop_back();
        
        LOG ("dfa-determinisation", 
             cout << "Pop worklist: " << ss_ID << " --> " << ss << "\n");
        
        for (unsigned int ia=alphastart; ia<alphastart+alphasz; ia++)
        {
          V termsymb = V::mkTerm(ia);
          set< unsigned int > dests;
          bool accept = false;
          for (set<unsigned int>::iterator ii=ss.begin(), ee=ss.end(); ii!=ee; ++ii)
          {
            this->getTransDest(*ii, termsymb, dests, accept);
          }
          
          LOG ("dfa-determinisation", 
               cout << "source: " << ss << " [" << termsymb << "] " 
               << " dest: " << dests;
               if (accept) cout << " (end state)\n"; else cout << endl);
          
          if (!dests.empty())
          {
            unsigned int nID = table.insert(dests, accept, ID);
            if (ID == nID) {
              ID++;	  
              worklist.push_back(nID);
              rev_table.insert(make_pair(nID,dests));
              
              LOG ("dfa-determinisation", 
                   cout << "Push worklist: " << nID << " --> " << dests << "\n");
              
            }
            dfa.transition(dfa.state(ss_ID), termsymb, dfa.state(nID));
            if (accept) dfa.accept(dfa.state(nID));
          }
          else
          { // to the sink state
            dfa.transition(dfa.state(ss_ID), termsymb, sink);
          }
        } // end outer for
      } // end while    
      
      // special case
      if (this->accepts[this->start]) 
        dfa.accepts[dfa.start] = true;
      
      // Transition for the sink state
      for (unsigned int ia=alphastart; ia<alphastart+alphasz; ia++)
        dfa.transition(sink, V::mkTerm(ia), sink);    
      
      dfa.eliminateUnreachableStates();
      
      LOG ("dfa-determinisation", 
           cout << "Conversion from NDFA to DFA done!\n");
      
      return dfa;
    }

    // Pre: *this is a deterministic finite automata
    // Swap initial and final states.
    void complement()
    {
      for(unsigned int i=0; i < this->accepts.size(); i++)
      {
        this->accepts[i].flip(); 
      }
    }

    // *this is the automata that recognizes the union of *this and
    // other.
    // Post: *this is a nondeterministic finite automata but
    //       without epsilon transitions.
    void nondet_union(const dfa_t &other)
    { 
      
      if (other.accepts.empty()) 
        return;
      if (other.trans.empty() && other.eps.empty()) 
        return;
      
      MapStates mapping; // mapping to relabel states old -> new
      
      unsigned int nstates = (unsigned int) this->nStates();
      for(int i=0; i< other.nStates(); i++)
      {
        mapping.insert(make_pair(i, nstates));
        this->state(nstates++);
      }
      
      for( unsigned int ss = 0; ss < other.trans.size(); ss++ )
      {
        for( unsigned int ti = 0; ti < other.trans[ss].size(); ti++ )
        {
          this->transition(mkState(mapping[ss]), 
                           other.trans[ss][ti].val,
                           mkState(mapping[other.trans[ss][ti].dest]));
        }
        for(unsigned int ei = 0; ei < other.eps[ss].size(); ei++)
        {
          this->eps_transition(mkState(mapping[ss]), 
                               mkState(mapping[other.eps[ss][ei]])); 
        }
      }
      
      // set new accepting states
      for (int i=0; i < other.nStates() ; i++)
      {
        if (other.accepts[i])
          this->accept(mkState(mapping[i]));
      }
      
      // finally, set new initial state
      State new_start       = this->state(this->nStates());
      State old_start_ndfa  = this->startState();
      State old_start_other = mkState(mapping[other.startState().id]);
      this->eps_transition(new_start, old_start_ndfa);
      this->eps_transition(new_start, old_start_other);
      this->setStart(new_start);
      // Make sure no epsilon transitions after union
      this->eps_elim();    
      bool is_not_empty = this->eliminateDeadStates();
      if (!is_not_empty)
        throw error("nondet union of DFA's returned an empty automata");
    }

    // Return the automata that recognizes the union of dfa1 and dfa2.
    // Pre : dfa1 and dfa2 are deterministic finite automata 
    // Post: the returned value is a deterministic finite automata.
    static dfa_t det_union(const dfa_t &dfa1, const dfa_t &dfa2)
    { 
      return product<UnionAccept>(dfa1,dfa2);
    }

    // Return true iff the intersection is not empty and inter is the
    // automata that recognizes the intersection between dfa1 and dfa2.
    // Pre: dfa1 and dfa2 are deterministic finite automata 
    static boost::optional<dfa_t> intersection(const dfa_t &dfa1, 
                                               const dfa_t &dfa2)
    { 
      dfa_t inter = product<IntersectAccept>(dfa1,dfa2);
      bool res = inter.eliminateDeadStates();
      if (res)
        return boost::optional<dfa_t> (inter);
      else
        return boost::optional<dfa_t> ();
    }
       
    // L(dif_sfa) = L(sfa1) \ L(sfa2).
    // Return true if  L(dif_sfa) is not empty.
    static boost::optional<dfa_t> difference(const dfa_t &fa1, 
                                             const dfa_t &fa2, 
                                             unsigned int alphstart, 
                                             unsigned int alphsz)
    {
      dfa_t dfa2 = fa2.makeDFA (alphstart, alphsz);
      dfa2.complement();
      dfa_t dfa1 = fa1.makeDFA(alphstart,alphsz);
      return intersection(dfa1, dfa2);
    }


    // An expensive but pretty simple minimization algorithm.
    // PRE : dfa must be a deterministic finite automata
    // POST: return the smallest equivalent automata
    dfa_t minimize(const unsigned int alphstart, const unsigned int alphsz)
    {

      dfa_t dfa(*this);

      // First eliminate unreachable states
      dfa.eliminateUnreachableStates();

      // Fill out the distinct table
      DFA_impl::distinctTable<bool> is_distinct(dfa.nStates(), false);
      for(int i=0; i < dfa.nStates(); i++)
      {
        for(int j=0; j < dfa.nStates(); j++)
        {
          if ( ( dfa.accepts[i] && !dfa.accepts[j]) ||
               (!dfa.accepts[i] &&  dfa.accepts[j]))
          {
            is_distinct[i][j] = true;
            //cout << "Marked " << "[" << i << "," << j << "]" << endl;
          }
        }
      } 
      bool changed = true;
      while (changed)
      {
        changed = false;
        for(int i=0; i < dfa.nStates(); i++)
        {
          for(int j=0; j < dfa.nStates(); j++)
          {
            if (!is_distinct[i][j]){
              bool marked=false;
              for (unsigned int k=alphstart; (k < alphstart+alphsz && !marked);k++)
              {
                if (is_distinct[dfa.delta (i,k)][dfa.delta (j,k)])
                {
                  marked = changed = is_distinct[i][j] = true;
                  //cout << "Marked transitively " << "[" << i << "," << j << "]" << endl;
                }
              }
            }
          }
        } 
      }

      ////
      // Create the minimal dfa
      ////
      
      dfa_t min_dfa;    
      // Create the equivalence classes of the "indistinguishability"
      // relation
      UnionFind eq_classes(dfa.nStates());
      for(int i=0; i < dfa.nStates(); i++)
      {
        for(int j=dfa.nStates()-1; j > i ; j--)
        {
          if (!is_distinct[i][j])
            eq_classes.set(i,j);
        }
      }
      // Create mapping to relabel states old->new
      MapStates mapping;
      unsigned new_state_id =0;
      vector<unsigned> new_states;
      for(int i=0; i < dfa.nStates();i++)
      {
        // we map only representatives 
        if (eq_classes.lookup(i) == i)
        {
          mapping.insert(make_pair(i,new_state_id++));
          new_states.push_back(i);
        }
      }
      
      // Create new states
      for(unsigned int i=0; i< new_states.size(); i++)
        min_dfa.state(i);
      // Marking initial state
      unsigned new_start = mapping[eq_classes.lookup(dfa.start)];
      min_dfa.setStart(mkState(new_start));
      // Copy delta transitions
      for (unsigned int i=0; i < new_states.size(); i++)
      {
        unsigned src = new_states[i];
        for(unsigned int j=0; j< dfa.trans[src].size(); j++)
        {      
          unsigned new_src = mapping[src];
          unsigned new_dst = mapping[eq_classes.lookup(dfa.trans[src][j].dest)];
          min_dfa.transition(mkState(new_src), 
                             dfa.trans[src][j].val, 
                             mkState(new_dst));
        }
      }
      // Marking final states
      for (unsigned int i=0; i < new_states.size(); i++)
      {
        unsigned s = new_states[i];
        if (dfa.accepts[s])
          min_dfa.accept(mkState(mapping[s]));        
      }
      return min_dfa;
    }
    
   private:

    // Pre: *this does not have epsilon transitions
    bool  eliminateUnreachableStates()
    {
      // initialize vector
      vector<bool> reachable_states;
      reachable_states.reserve(this->nStates());
      for(int i=0; i< this->nStates(); i++)
        reachable_states.push_back(false);
      
      // mark every state reachable from the start state.
      vector<unsigned int> worklist;
      worklist.push_back(this->startState().id);
      reachable_states[this->startState().id] = true;
      while (!worklist.empty())
      {
        unsigned int qi = worklist.back();
        worklist.pop_back();
        for (unsigned int qj=0; qj < this->trans[qi].size(); qj++)
        {
          if (!reachable_states[this->trans[qi][qj].dest])
          {
            reachable_states[this->trans[qi][qj].dest] = true;
            worklist.push_back(this->trans[qi][qj].dest);
          }
        }
      } 
      bool is_accepting;
      *this = resize(reachable_states, is_accepting);
      return is_accepting;
    }

    // Pre: *this does not have epsilon transitions
    // Pre: the size of states must be the number of states in *this.
    // Return true iff there is a path from qi to some accepting state.
    // It performs a simple depth-first search.
    bool mark_generating_states(State qi, 
                                set<State> &visited, 
                                vector<bool> &states) const
    {
      if (visited.count(qi)) 
        return states[qi.id];

      visited.insert(qi);

      if (this->trans[qi.id].empty())
      {
        states[qi.id] = this->accepts[qi.id];
      }
      else
      {
        bool is_generating=false;
        for (unsigned int qj=0; qj < this->trans[qi.id].size(); qj++)
        {
          is_generating |= mark_generating_states(mkState(this->trans[qi.id][qj].dest), 
                                                  visited, 
                                                  states);
        }
        states[qi.id]= (this->accepts[qi.id] | is_generating);
      }
      return states[qi.id];
    }
    
    // Pre: *this has no epsilon transitions
    // Return true iff there is a path from qi to any accepting state.
    bool is_generating_state(State qi) const
    { 
      set<State> visited;
      vector<bool>    states;
      states.reserve(this->nStates());
      for(int i=0; i< this->nStates(); i++)
      {
        states.push_back(false);
      }
      return mark_generating_states(qi, visited, states);
    }
            
    // Helper for makeDFA
    void getTransDest(unsigned int src, 
                      V symbol, 
                      set<unsigned int> &dests, 
                      bool &accept) const
    {
      for (unsigned int i=0; i < this->trans[src].size() ; i++)
      {
        if (this->trans[src][i].val == symbol)
        {
          dests.insert(this->trans[src][i].dest);
          accept |= this->accepts[this->trans[src][i].dest];
        }
      }
    }

    // Helper for product construction (see below)
    // Add a new state <s1,s2> in sfa and return the new state.
    State mkPairState(unsigned int s1, unsigned int s2, 
                      bool accept,
                      DFA_impl::StateTable<vector<unsigned int> > &table, 
                      unsigned int &ID)
    {
      
      vector<unsigned int> newState;
      newState.push_back(s1); newState.push_back(s2);
      unsigned int nID = table.insert(newState, accept , ID);
      if (ID == nID)
      {  // new state
        
        LOG ("dfa-product", 
             cout << "Inserted new " << newState << " final:" << accept
             << " --> " <<  nID << endl);
        
        State q = this->state(ID);
        ID++;
        if (accept) 
        {
          LOG ("dfa-product", 
               cout << "Marked " << nID << " as final state\n");
          
          this->accept(q);
        }
        return q;
      }
      else
      {						
        
        LOG ("dfa-product", 
             cout << "Reused " << newState << " final:" << accept
             << " --> " <<  nID << endl);
        
        return mkState(nID);
      }
    }

    // Construct the product automata using classical product construction
    // Pre: sfa1 and sfa2 do not have epsilon transitions
    // Post: product is the automata that recognizes the intersection or
    // union (depends on the template parameter T) of the languages of
    // sfa1 and sfa2. 
    template<typename T>
    static dfa_t product(const dfa_t &sfa1, const dfa_t &sfa2)
    { 
      LOG ("dfa-product", 
           cout << "Begin product construction of sfa's ... \n");
      
#ifdef SANITY_CHECKS    
      sfa1.CheckNoEpsilonTrans();
      sfa2.CheckNoEpsilonTrans();
#endif 

      dfa_t product;

      // Initialise bookkeeping.
      unsigned int ID = 0;
      DFA_impl::StateTable<vector<unsigned int> >  table;
      table.clear();
      // Create all states from all pairs 
      // Also, marking final and initial states.
      for (unsigned int q1 = 0; q1 < (unsigned int) sfa1.nStates(); q1++){
        for (unsigned int q2 = 0; q2 < (unsigned int) sfa2.nStates(); q2++){
          T acceptor;
          const bool accept = acceptor(sfa1.accepts[q1], sfa2.accepts[q2]);
          State q = product.mkPairState(q1, q2, accept , table,  ID);
          if (q1 == sfa1.startState().id && q2 == sfa2.startState().id){
            
            LOG ("dfa-product", 
                 cout << "Marked " << q.id << " as initial state\n");
            
            product.setStart(q);
          }
        }
      }
      
      LOG ("dfa-product",   cout << endl);
      
      // Create delta transitions
      for (unsigned int q1_src = 0; q1_src < (unsigned int) sfa1.nStates(); q1_src++)
      {
        for (unsigned int q2_src = 0; q2_src < (unsigned int) sfa2.nStates(); q2_src++)
        {
          for(unsigned int ei = 0 ; ei < sfa1.trans[q1_src].size() ; ei++)
          {
            for(unsigned int ej = 0 ; ej < sfa2.trans[q2_src].size() ; ej++)
            {
              if (sfa1.trans[q1_src][ei].val == sfa2.trans[q2_src][ej].val)
              {
                unsigned int q1_dest = sfa1.trans[q1_src][ei].dest ;
                unsigned int q2_dest = sfa2.trans[q2_src][ej].dest ;
                
                LOG ("dfa-product",
                     cout << "SFA1: "
                     << q1_src << " ---> " << q1_dest << " with " 
                     << sfa1.trans[q1_src][ei].val << endl;
                     cout << "SFA2: "
                     << q2_src << " ---> " << q2_dest << " with " 
                     << sfa2.trans[q2_src][ej].val << endl);
                
                T acceptor;
                const bool accept_src  = acceptor (sfa1.accepts[q1_src], 
                                                   sfa2.accepts[q2_src]);
                const bool accept_dest = acceptor (sfa1.accepts[q1_dest], 
                                                   sfa2.accepts[q2_dest]);

                State src = product.mkPairState(q1_src, q2_src, accept_src , table, ID);
                State dest= product.mkPairState(q1_dest, q2_dest, accept_dest , table, ID);
                                    
                LOG ("dfa-product",
                     cout << "PRODUCT: " << src.id << "=<" << q1_src << "," << q2_src << ">" 
                     << " ---> " << dest.id << "=<" << q1_dest << "," << q2_dest << ">" 
                     <<  " with " << sfa1.trans[q1_src][ei].val << endl);
                
                product.transition(src, sfa1.trans[q1_src][ei].val, dest);
                
                LOG ("dfa-product", cout << endl);
              }
            }
          }
        }
      } 
      product.eliminateUnreachableStates();
      
      LOG ("dfa-product",
           cout << "End product construction of sfa's. \n" << product << "\n");

      return product;
    }

    struct IntersectAccept
    {
      bool operator()(bool a, bool b){ return a & b;}
    };
    
    struct UnionAccept
    {
      bool operator()(bool a, bool b){ return a | b;}
    };


    // Helper to apply the delta function for minimize
    // It's not constant!
    unsigned delta(const unsigned p, const unsigned sym) const
    {
      for( unsigned int ti = 0; ti < this->trans[p].size(); ti++ )
      {
        if (this->trans[p][ti].val.symID () ==  (int) sym)
          return this->trans[p][ti].dest;
      }
      throw error("DFA is not well formed");
    } 

    
  }; // end class DFA

  template<class V>
  inline std::ostream& operator<<(std::ostream& o, DFA<V> dfa) 
  {
    return dfa.write(o);
  }

  // Calculates the shortest path from each state to an accept, then
  // returns the longest.
  template<class V>
  int max_minpath(DFA<V>& dfa)
  {
    // The set of inverted edges.
    std::vector< std::vector<int> > inverse(dfa.nStates());

    // Breadth-first queue of states to be processed
    SparseSet<> queue(dfa.nStates());
    // Distance from state i to an accept
    std::vector<int> dist(dfa.nStates(), -1);

    for(unsigned int si = 0; si < dfa.trans.size(); si++)
    {
      if(dfa.accepts[si])
      {
        queue.insert(si);
        dist[si] = 0;
      }

      for(unsigned int ei = 0; ei < dfa.trans[si].size(); ei++)
      {
        inverse[dfa.trans[si][ei].dest].push_back(si);
      }
    }

    unsigned int qi = 0;
    for(; qi < queue.size(); qi++)
    {
      int si = queue[qi];
      int sd = dist[si];
      std::vector<int> state_inv(inverse[si]);
      for(unsigned int ei = 0; ei < state_inv.size(); ei++)
      {
        if(!queue.elem(state_inv[ei]))
        {
          queue.insert(state_inv[ei]);
          dist[state_inv[ei]] = sd+1;
        }
      }
    }
    return dist[queue[qi-1]];
  }

} // end namespace covenant
#endif
