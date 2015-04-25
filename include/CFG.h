#ifndef __CFG_H__
#define __CFG_H__

#include <iostream>
#include <stack>
#include <limits>

#include <boost/unordered_map.hpp>

#include <Common.h>
#include <adt/MapSet.h>

using namespace std;

namespace covenant 
{

  class Sym 
  {
    int x;
    
   public:
    
    static Sym mkTerm(int id)
    {
      Sym p;
      p.x = (id<<1); 
      return p;
    }

    static Sym mkVar(int id)
    {
      Sym p; 
      p.x = (id<<1) + 1; 
      return p;
    }
    
    bool operator == (Sym p) const { return x == p.x; }
    bool operator != (Sym p) const { return x != p.x; }
    bool operator < (Sym p)  const { return x < p.x; }

    bool isVar()  const  { return this->x&1; }
    bool isTerm() const  { return !this->isVar(); }
    int  symID()  const  { return this->x>>1; }
    
    friend ostream& operator<<(ostream& o, Sym s) 
    {
      if( s.isVar () )
        s.print_var(o);
      else 
        o << s.symID ();
      return o;
    }

   private:

    void print_var(ostream& o)
    {
      o << "NT_" << symID();
    }
  };  

  
  struct HashSym 
  {
    size_t operator()(const Sym &s) const
    {
      boost::hash<const int> hasher;
      return hasher(s.symID());
    }
  };

  struct EqSym 
  {
    bool operator()(const Sym &x, const Sym &y) const
    {
      return x.symID () == y.symID ();
    }
  };

  
  inline bool IsTerm(Sym s)
  {
    return s.isTerm ();
  }

  // Succeed if all symbols in S are terminal 
  inline bool isAllTerminal(const vector<Sym> &S) 
  {
    unsigned numTerm = count_if(S.begin(),S.end(), IsTerm);
    return (numTerm == S.size());
  }

  class Rule 
  {

    TermFactory _tfac;

   public:

    vector<Sym> _syms;

    Rule (TermFactory tfac): 
        _tfac (tfac) { }

    explicit Rule (const vector<Sym> &syms, TermFactory tfac): 
        _tfac (tfac), _syms(syms) { }

    Rule (const Rule &other): 
      _tfac(other._tfac), _syms(other._syms)  { }
  
    Rule& operator=(Rule other)
    {
      this->_tfac = other._tfac;
      this->_syms = other._syms;
      return *this;
    }

    // Make an epsilon rule
    static Rule E (TermFactory tfac) 
    { 
      return Rule (tfac); 
    }

    TermFactory getTermFactory() {  return _tfac; }

    Rule& operator<< (Sym s)
    {
      _syms.push_back(s);
      return *this;
    }
    
    Rule& operator<< (int n)
    {
      _syms.push_back(Sym::mkTerm(n));
      return *this;
    }
    
    Sym& operator[] (unsigned int i) { return _syms[i]; }

    bool operator== (const Rule &other) const 
    { 
      if (other._syms.size() != this->_syms.size()) 
        return false;

      return equal(other._syms.begin(), 
                   other._syms.end(), 
                   this->_syms.begin());
    }

    friend ostream& operator<<(ostream& o, Rule r)
    {
      if( r._syms.size() == 0 )
      { 
        o << "e"; 
      } 
      else 
      {
        TermFactory tfac = r.getTermFactory ();
        assert (tfac);
        for( unsigned int i=0; i < r._syms.size(); i++ )
        {
          if (r._syms[i].isVar ())
          {
            o << r._syms[i] << " " ;
          }
          else
          {
            // we apply the dictionary to output the corresponding
            // string.
            o << tfac->remap(r._syms[i].symID ()) << " " ;
          }
        }
      }
      return o;
    }
    
  }; 


  // FIXME: 
  // - make Sym a template parameter
  // - replace the use of vector<Sym> with Rule
  // - make private class attributes 
  //  - use iterators to traverse grammar productions
  class CFG 
  {

    typedef boost::unordered_map<unsigned int, unsigned int> term_rule_map_t;
    term_rule_map_t TrackedTerms;

  public:

    struct ProdInfo 
    {
      ProdInfo(int _rule): rule(_rule) { } 
      int rule;
    }; /* end class ProdInfo */

    int alphsz;    // size of the CFG alphabet
    int alphstart; // we assume that the CFG alphabet ranges from {alphstart, ..., alphsz-1}
    int start;     // Start symbol of the CFG

    vector< vector<ProdInfo> > prods; // Mapping var -> rule.
    vector< vector<Sym> > rules;      // The actual production rules.

    TermFactory _tfac;

  public:

    CFG (TermFactory tfac): 
        alphsz(std::numeric_limits<short int>::max ()), 
        alphstart(0), 
        start(0),
       _tfac(tfac)
    { }

    CFG (int _alphsz, TermFactory tfac): 
        alphsz(_alphsz), 
        alphstart(0), 
        start(0),
        _tfac(tfac) 
    { }

    CFG (const CFG& other): 
      TrackedTerms (other.TrackedTerms),
      alphsz(other.alphsz), alphstart(other.alphstart), start(other.start),
      prods (other.prods), rules (other.rules), _tfac(other._tfac) { }

    CFG& operator=(CFG other)
    {
      if (this != &other)
      {
        TrackedTerms  = other.TrackedTerms;
        alphsz = other.alphsz;
        alphstart = other.alphstart;
        start = other.start;
        prods  = other.prods;
        rules = other.rules;
        _tfac = other._tfac;
      }
      return *this;
    }
                           

    ~CFG () 
    {
      prods.clear ();
      rules.clear ();
      TrackedTerms.clear ();
    }
    
    TermFactory getTermFactory () { return _tfac; }

    Sym term(int i) 
    { 
      return Sym::mkTerm(i); 
    }

    Sym newVar()
    {
      int id = prods.size();
      prods.push_back( vector<ProdInfo>());
      return Sym::mkVar (id);
    }
    
    void setStart(Sym s)
    {
      if (!s.isVar()) 
        throw error("start symbol of the CFG must be nonterminal");

      start = s.symID(); 
    }
    
    Sym startSym(void) const 
    {
      return Sym::mkVar(start);
    }
    
    int nVars(void) const 
    {
      return prods.size();
    }
    
    // Make a grammar production
    void prod(Sym LHS, const Rule& RHS)
    {
      int r_id = rules.size();
      rules.push_back(RHS._syms); // Add the grammar rule
      prods[LHS.symID ()].push_back(ProdInfo(r_id)); // Update the mapping
    }
    
    // Transform production rules to ensure that:
    // - there are at most two symbols in each rule,
    // - the lhs of a production does not appear on the rhs, and
    // - if |rhs| > 1 then it does not contain any terminal symbol.
    void normalize(void)
    {
      for ( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
	for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
          {
	  const int r = prods[vv][ri].rule;

            // Step 1: at most two symbols in each production
            vector<Sym> nr(2);
            while( rules[r].size() > 2 )
            {
              Sym next( newVar() );
              nr[1] = replaceTerminal(rules[r].back());
              rules[r].pop_back();
              nr[0] = replaceTerminal(rules[r].back());
              rules[r].pop_back();
              rules[r].push_back(next);
              rules.push_back(nr);
              prods[next.symID ()].push_back(ProdInfo(rules.size()-1));
            }

            // Step 2: lhs of the production does not appear on the rhs
            vector<Sym> SymVec(rules[r]);
            replaceNonTerminal (SymVec, Sym::mkVar(vv));
            rules[r] = SymVec;

            assert (rules[r].size() <= 2);
            
	  // Step 3: production of length 2 cannot have nonterminals
	  if (rules[r].size() == 2)
            {
              vector<Sym> SymVec(rules[r]);
              replaceTwoUnitWithMemoing (SymVec);
              rules [r] = SymVec;
            }

	} // end inner for 
      } // end outer for      

    }

    // PRE: the grammar is normalized
    void get_alphabet_min_and_max(int &min, int &max)
    {
      vector<int> terminals;
      for ( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
        for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
        {
          const int r = prods[vv][ri].rule;
          vector<Sym> RHS = rules[r];
          if (RHS.size() == 1 && RHS[0].isTerm())
            terminals.push_back(RHS[0].symID());
        }
      }
      min = *min_element(terminals.begin(), terminals.end());
      max = *max_element(terminals.begin(), terminals.end());
    }

    void eliminateEpsilon()
    {
      // First pass to add into NullableSymbols initial epsilon rules.
      // NullableSymbols contains the lhs of A -> epsilon
      set<Sym> NullableSymbols;
      for ( unsigned int vv = 0; vv < prods.size(); vv++ ){
	for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ ){
	  const int r = prods[vv][ri].rule;
	  if (rules[r].size() == 0) { // epsilon rule
	    Sym nt = Sym::mkVar(vv);
	    NullableSymbols.insert(nt);
	    // delete already the epsilon rule
	    // cout << "Deleting epsilon rule " << nt << " -> " << rules[r] << "\n";
	    removeRule(vv,ri);
	  }
	} // end for
      } //end for
      
      // Second pass: not sure if we really need two passes.

      // FIXME: add a while until no change (new epsilon rules)
      // Also, if we have only rule A -> epsilon we can remove A on all rhs directly!
      // (should be done in replaceNullableSymol)
      vector<pair<Sym,Rule> > newRules;
      for( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
	for( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
          {
	  const int r = prods[vv][ri].rule;
	  vector<Sym> OldRHS = rules[r];
	  // Add indirect epsilon rules
	  if ((OldRHS.size() == 1) &&  
                (NullableSymbols.count(OldRHS[0]) > 0))
            {
	    removeRule(vv,ri);
	    NullableSymbols.insert(Sym::mkVar(vv));
	  }
	  else{  
	    set<vector<Sym> > tmpRHSs = 
	      replaceNullableSymbol(OldRHS, NullableSymbols);
	    typedef set<vector<Sym> >::iterator It;
	    for (It it1 = tmpRHSs.begin(), et1 = tmpRHSs.end() ; it1 != et1; ++it1)
              {
	      //Add new production: do not add directly into the grammar
	      //or infinite loop.
	      newRules.push_back(make_pair(Sym::mkVar(vv),
                                             Rule(*it1, _tfac)));
	    }
	  }
	}
      }
      // Add into the grammar the new rules
      while (!newRules.empty())
      {
        pair<Sym,Rule> r =  newRules.back();
        if (!isRedundantRule(r.first,r.second))
        {
          prod(r.first,r.second);
        }
        newRules.pop_back();
      }
    }

    // Return true if the grammar rules are in Chomsky Normal Form. 
    // That is, if all rules are of the form either A -> BC, or S ->
    // epsilon, or A -> alpha.
    bool IsCNF() const 
    {
      for ( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
        for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
        {
          const int r = prods[vv][ri].rule;
          vector<Sym> RHS = rules[r];
          // We only allow S -> epsilon where S is start symbol
          if( RHS.size() == 0 && (startSym() != Sym::mkVar(vv)))
            return false; 
          if (RHS.size() == 1 && RHS[0].isVar ())
            return false;
          if (RHS.size() == 2 && !(RHS[0].isVar () && RHS[1].isVar ()))
            return false;
          if (RHS.size() > 2)
            return false;
        }
      }
      return true;
    }

    inline bool isRegularGrammar() const
    {
      // Note that we cannot mix left and right regular rules.
      // E.g., G = S->aA; A->Sb; S->e contains only either left or
      // right regular rules. However, L(G) is not regular!
      return (isLeftRegGrammar() || isRightRegGrammar());
    }

    // Convert to Chomsky normal form. 
    // Production can be only of the form:
    // S -> epsilon
    // A -> BC | alpha
    void ConvertToCNF()
    {
      // new start symbol
      Sym newStart = newVar();
      Sym oldStart = startSym();
      prod(newStart, Rule::E(_tfac) << oldStart);
      setStart(newStart);      

      LOG ("cfg-cnf", cout << "After adding new start: " << endl << *this);

      normalize();

      LOG ("cfg-cnf", cout << "After normalization: " << endl << *this);

      eliminateEpsilon();

      LOG ("cfg-cnf", cout << "After epsilon elimination: " << endl << *this);

      eliminateUnit();
      
      LOG ("cfg-cnf", cout << "After unit elimination: " << endl << *this);

      if (!IsCNF())
      {
        throw error("CFG is not in CNF");
      }
      // TODO: eliminate useless symbols: non-generating symbols, and
      // unreachable symbols.
    }
    
    set<Sym> computeGenerating()
    {
      set<Sym> GenSet;
      bool change=true;
      while (change)
      {
        change = false;
        for ( unsigned int vv = 0; vv < prods.size(); vv++ )
        {
          for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
          {
            const int r = prods[vv][ri].rule;
            vector<Sym> RHS = rules[r];
            bool allSatisfy=true;
            for(unsigned si =0; (si < RHS.size() && allSatisfy) ; si++)
            {
              if (!(RHS[si].isTerm() || GenSet.count(RHS[si]) > 0))
                allSatisfy=false;
            }
            if (allSatisfy)
              change |= InsertAndNotifyChange(GenSet,Sym::mkVar(vv));
          }
        }
      }
      return GenSet;
    }

    // Remove first all non-generating variables along with
    // productions that involve non-generating variables, and second
    // remove unreachable variables along with productions that
    // involve unreachable variables.
    // E.g., given the CFG:
    // S -> AB | a
    // A -> aA
    // B -> b
    // and given GenSet={S,B}, ReachSet={S,A,B} we produce a new CFG:
    // S -> a   
    void minimize()
    {
      throw error("minimization of CFGs is not implemented");
    }

    inline void check_well_formed()
    {
      bool error=false;
      ostringstream error_msg;
      if (prods.size() == 0)
      {
        error=true;
        error_msg << "ERROR: grammar must have at least one production rule." << endl;
      }
      for( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
        if (prods[vv].size() == 0){
          error=true;
          error_msg << "ERROR: " 
                    << Sym::mkVar(vv) << " non-terminal symbol without rule.\n"; 
          error_msg << "Also possible a terminal misinterpreted as nonterminal.\n";
        }
      }
      for ( unsigned int vv = 0; (vv < prods.size() && !error); vv++ )
      {
        for ( unsigned int ri = 0; (ri < prods[vv].size() && !error); ri++ )
        {
          const int r = prods[vv][ri].rule;
          vector<Sym> RHS = rules[r];
          if (RHS.size() == 0) continue;
          for (unsigned k=0; ( k < RHS.size() && !error) ;k++)
          {
            if (RHS[k].isTerm () && 
                !(RHS[k].symID () >= alphstart && 
                  RHS[k].symID()  <  alphstart+alphsz)){
              error=true;
              error_msg << "ERROR: some terminal symbol is out of range. ";
              error_msg << "It must be between " << alphstart << " and " << alphstart+alphsz << endl;   
            }
          }
        }
      }
      if (error)
      {
        cout << *this << endl << error_msg.str();
        exit(1);
      }
    }

    ostream & write(ostream &o) const 
    {
      o << "start symbol: " << startSym () << endl;
      for( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
        printProduction(o, vv);
      }
      return o;
    }

    void get_stats(unsigned& terminals,
                   unsigned& nonterminals,
                   unsigned& productions) const
    {
      terminals = alphsz;
      nonterminals = prods.size ();
      productions = 0;
      for ( unsigned vv = 0; vv < prods.size(); ++vv)
        productions += prods[vv].size();
    }


    void stats(ostream &o) const
    {
      o << "Number of terminals        : " << alphsz << "\n";
      o << "Number of non-terminals    : " << prods.size () << "\n";

      unsigned total_productions=0;
      unsigned max_per_nonterminal=0;
      for ( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
        total_productions += prods[vv].size();
        max_per_nonterminal = std::max(max_per_nonterminal, 
                                       (unsigned) prods[vv].size ());
      }

      o << "Total number of productions: " << total_productions << "\n";
      if (prods.size() > 0)
      {
        o << "Avg number of productions per non-terminal : " 
          << float(total_productions) / float(prods.size ()) << "\n";
        o << "Max number of productions per non-terminal : " 
          << max_per_nonterminal << "\n";
      }
    }

    friend ostream& operator<<(ostream &o, CFG cfg)
    {
      cfg.write(o);
      return o;
    }


   private:

    // Replace a terminal with a fresh nonterminal symbol
    Sym replaceTerminal (Sym t)
    {
      if (t.isVar ()) return t;

      Sym nt = newVar();
      prod (nt, Rule::E(_tfac) << t);
      return nt;
    }

    // Replace a terminal with a nonterminal but trying to reuse an
    // existing nonterminal before creating a fresh one.
    Sym replaceTerminalWithMemoing(Sym t)
    {
      if (t.isVar ()) return t;
    
      term_rule_map_t::iterator It = TrackedTerms.find(t.symID ());
      if (It != TrackedTerms.end())
        return Sym::mkVar((*It).second);
      else
      {
        Sym nt = replaceTerminal(t);
        TrackedTerms.insert(make_pair(t.symID (), nt.symID ()));
        return nt;
      }
    }

    // Replace a terminal with a fresh nonterminal symbol
    void replaceTwoUnitWithMemoing(vector<Sym> & SymVec)
    {
      assert(SymVec.size () == 2);
      Sym a(SymVec[0]);
      Sym b(SymVec[1]);

      SymVec.clear ();
      SymVec.push_back (replaceTerminalWithMemoing (a));
      SymVec.push_back (replaceTerminalWithMemoing (b));
    }

    // replace a nonterminal with a fresh nonterminal
    void replaceNonTerminal (vector<Sym> &SymVec, Sym s)
    {
      for (unsigned int i=0; i<SymVec.size (); ++i)
      {
        if (SymVec [i].isVar () && (SymVec [i] == s))
        {
          Sym FreshNonTerm = newVar();
          prod (FreshNonTerm, Rule::E (_tfac) << SymVec [i]);
          SymVec [i] = FreshNonTerm;
        }
      }
    }

    vector<unsigned>  
    getNullablePositions(const vector<Sym> &SymVec, 
                         const set<Sym> &NullableSymbols)
    {
      vector<unsigned> positions;
      for (unsigned i=0; i < SymVec.size() ; i++)
      {		
        if (NullableSymbols.count(SymVec[i]) > 0)
          positions.push_back(i);
      }
      return positions;
    }

    set<vector<Sym> > 
    replaceNullableSymbol(const vector<Sym> &SymVec,
                          const set<Sym> &NullableSymbols)
    {
      typedef pair <vector<Sym>, vector<unsigned> > SymPosPair;
      
      set<vector<Sym> > newRules;
      vector<unsigned> positions = getNullablePositions(SymVec, NullableSymbols);
	
      if (positions.empty()) // no symbols with epsilon rules
	return newRules;

      vector<SymPosPair> worklist;
      worklist.push_back(make_pair(SymVec,positions));
      // BCB where B->epsilon produces {BCB, BC, C, CB}
      while (!worklist.empty())
      {
        SymPosPair x = worklist.back();
        worklist.pop_back();
        // cout << x.first << " " << x.second  << endl;
        if (x.second.empty()) // no more positions to replace
          newRules.insert(x.first);
        else
        {
          unsigned pos = x.second.back(); 
          x.second.pop_back(); 
          worklist.push_back(x); // added without replacement
          vector<Sym> new_rhs(x.first);
          remove(new_rhs,pos);   // added with replacement
          worklist.push_back(make_pair(new_rhs,
                                       getNullablePositions(new_rhs,
                                                            NullableSymbols)));
        }
      }
      return newRules;
    }

    // To avoid duplicated rules
    bool isRedundantRule(Sym lhs, const Rule &rhs) const
    {
      unsigned vv = lhs.symID ();
      for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
      {
        const int r = prods[vv][ri].rule;
        if (rhs == Rule (rules[r], _tfac))
          return true;
      }
      return false;
    }

    // prods[0][0] = 0
    // prods[0][1] = 1
    // prods[1][0] = 2 <----
    // prods[2][0] = 3
    // rules[0] = AB
    // rules[1] = a
    // rules[2] = aA
    // rules[3] = b
    void removeRule(unsigned i, unsigned j )
    {
      const int r = prods[i][j].rule;

      remove(prods[i],j);
      remove(rules, r);

      // update the mappings in prods
      for ( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
        for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
        {
          if (prods[vv][ri].rule > r)
            prods[vv][ri].rule--;
        }
      }
    }

    vector<Rule> getRules(Sym lhs)
    {
      unsigned vv = lhs.symID ();
      vector<Rule> new_rules;
      for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
      {
        const int r = prods[vv][ri].rule;
        new_rules.push_back(Rule(rules[r], _tfac));
      }
      return new_rules;
    }

    void eliminateUnit()
    {
      vector<pair<unsigned,vector<Sym> > > UnitRules;
      for ( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
	for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
          {
	  const int r = prods[vv][ri].rule;
	  vector<Sym> unit_rule = rules[r];
	  if ((unit_rule.size() == 1) && unit_rule[0].isVar ())
	    UnitRules.push_back(make_pair(vv,unit_rule));
	}
      }
      // A -> BCB|C|CB|BC
      // B -> 12
      // C -> A
      // Unit rules: 
      // A -> C
      // C -> A
      // to eliminate A->C whenever C -> w appears add rule A -> w
      // A -> BCB|CB|BC|A
      // B -> 12

      set<pair<unsigned,vector<Sym> > > ProcessedUnitRules;

      while (!UnitRules.empty())
      {
        pair<unsigned,vector<Sym> > U = UnitRules.back();
        UnitRules.pop_back();
        unsigned vv = U.first;
        vector<Sym> rhs = U.second;
        
        if (ProcessedUnitRules.count(make_pair(vv,rhs)) > 0)
          continue;
        
        Sym lhs = Sym::mkVar(vv);
        // eliminate unit rule
        
        // cout << "Eliminating unit rule: " << lhs << " -> " << rhs << endl;
        
        // for unit rule A->B whenever B -> w appears add rule A -> w
        vector<Rule> new_rules = getRules(rhs[0]);
        ProcessedUnitRules.insert(make_pair(vv,rhs));
        // add new rule and update worklist
        for(unsigned k=0; k < new_rules.size(); k++)
        {
          Rule r_k = new_rules[k];
          //  cout << "Adding new rule: " << lhs << " -> " << r_k._syms << endl;
          if (!isRedundantRule(lhs,r_k))
          {
            prod(lhs, r_k);  
            if ( (r_k._syms.size() == 1) && (r_k._syms[0].isVar ()))
              UnitRules.push_back(make_pair(vv,r_k._syms ));
          }
        }
      }
      // cout << "Before removing unit rules: " << endl;
      // cout << *this;
      // Finally, removing unit rules
      for ( unsigned int vv = 0; vv < prods.size(); vv++ )
      {
        for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
        {
          const int r = prods[vv][ri].rule;
          vector<Sym> unit_rule = rules[r];
          if ((unit_rule.size() == 1) && unit_rule[0].isVar ())
          {
            // cout << "removing " << vv << "-" << ri << endl;
            removeRule(vv,ri);
            // cout << *this;
          }
        }
      }
    }
    
    bool isOneLevelIndirectionTerminal(const Sym &s) const
    {
      if (s.isTerm ()) 
        return true;
      else
      {
        bool is_term = true;
        int vv = s.symID();
        for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
        {
          const int r = prods[vv][ri].rule;
          vector<Sym> RHS = rules[r];
          is_term &= isAllTerminal(RHS); 
        }
        return is_term;
      }
    }

    // PRE: the grammar must be first normalized
    // Return true if all grammar rules are of the form 
    //    A -> e | aB | a | BC and B -> a
    // where A,B,C are any non-terminals and  "a" is any terminal. 
    bool isLeftRegGrammar() const
    {
      bool is_left = true;
      for ( unsigned int vv = 0; (vv < prods.size() && is_left); vv++ )
      {
        for ( unsigned int ri = 0; (ri < prods[vv].size() && is_left); ri++ )
        {
          const int r = prods[vv][ri].rule;
          vector<Sym> RHS = rules[r];
          unsigned int size = RHS.size();
          if (size == 0) 
            continue;
          else if (size == 1)
            is_left &= isOneLevelIndirectionTerminal(RHS[0]);
          else if (size == 2)
            is_left &= (isOneLevelIndirectionTerminal(RHS[0]) && RHS[1].isVar());
          else 
            is_left = false;
        }
      }
      return is_left;
    }

    // PRE: the grammar must be first normalized
    // Return true if all grammar rules are of the form 
    //    A -> e | Ba | a | BC and B -> a
    // where A,B,C are any non-terminals and  "a" is any terminal. 
    bool isRightRegGrammar() const 
    {
      bool is_right = true;
      for ( unsigned int vv = 0; (vv < prods.size() && is_right); vv++ )
      {
        for ( unsigned int ri = 0; (ri < prods[vv].size() && is_right); ri++ )
        {
          const int r = prods[vv][ri].rule;
          vector<Sym> RHS = rules[r];
          unsigned int size = RHS.size();
          if (size == 0) 
            continue;
          else if (size == 1)
            is_right &= isOneLevelIndirectionTerminal(RHS[0]);
          else if (size == 2)
            is_right &= (isOneLevelIndirectionTerminal(RHS[1]) && RHS[0].isVar ());
          else 
            is_right = false;
        }
      }
      return is_right;
    }


    inline bool InsertAndNotifyChange(set<Sym> &S, Sym e)
    {
      if (S.count(e) > 0) { return false; }
      else{ S.insert(e); return true; }
    }

    set<Sym> computeReachable()
    {
      set<Sym> ReachSet;
      ReachSet.insert(startSym());
      bool change=true;
      while (change)
      {
        change=false;
        for ( unsigned int vv = 0; vv < prods.size(); vv++ )
        {
          for ( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
          {
            const int r = prods[vv][ri].rule;
            vector<Sym> RHS = rules[r];
            if (ReachSet.count(Sym::mkVar(vv)) > 0){
              for(unsigned si =0; si < RHS.size(); si++)
              {
                if (RHS[si].isVar ())
                  change |= InsertAndNotifyChange(ReachSet,RHS[si]);
              }
            }
          }
        }
      }
      return ReachSet;
    }

    void printProduction(ostream &o, unsigned vv) const 
    {
      o << Sym::mkVar(vv);
      o << " -> ";
      for( unsigned int ri = 0; ri < prods[vv].size(); ri++ )
      {
        if( ri != 0 ) 
        { 
          o << "|"; 
        }
        o << Rule(rules[prods[vv][ri].rule], _tfac );
      }
      cout << endl;
    }
    
  }; /* end CFG class */

  // Tarjan's algorithm implemented based on code from:
  // https://github.com/leegoex/CPPStrongConnection
  // which is essentially a verbatim version of the Wikipedia pseudocode.
  class CFGConnect  {
    typedef int VID;

    class VertInfo {
    public:
    VertInfo(void)
      : index(-1), lowlink(-1), group_id(-1)
	{ }
      int index;
      int lowlink;
      int group_id;
      vector<int> edges;
    };
    
    int index;
    stack<int> S;
    
    vector<VertInfo> verts;
    vector< vector<int> > groups;

  public:

    CFGConnect(): index(-1) { }

    ~CFGConnect(){ cleanup(); }
    
    void cleanup(){
      verts.clear();
      groups.clear();
      while (!S.empty()) { S.pop(); }
    }

    void scc(CFG& g){
      cleanup();
      SCC_Init(g);
      SCC_Compute(verts.size()-1);
      
      // Eliminate the fake vertex.
      verts.pop_back();
      groups.pop_back();

      // Sort to make queries more efficient
      for (unsigned int i=0; i < groups.size(); i++){
	sort(groups[i].begin(), groups[i].end());
      }
    }
    
    const unsigned int nGroups(void) const { return groups.size(); }

    const vector<int>& group(int g){
      assert(g < (int) groups.size());
      return groups[g];
    }
    
    int group_id(int v){
      return verts[v].group_id;
    }

    ostream& write(ostream& o) const{
      for(unsigned int ii = 0; ii < nGroups(); ii++){
        const vector<int>& group(groups[ii]);
        o << "[";
        if(group.size() > 0)
          o <<  "NT_" << group[0]; //((char) ('A' + group[0]));
        for(unsigned si = 1; si < group.size(); si++)
          o << ", " << "NT_" << group[si]; //((char) ('A' + group[si]));
        o << "] ";
      }
      return o;
    }

    friend ostream& operator<<(ostream& o, CFGConnect conn){
      conn.write(o);
      return o;
    }

  private:
    void SCC_Compute(VID vi){
      VertInfo& v(verts[vi]);
      v.index = index;
      v.lowlink = index;
      index += 1;
      S.push(vi);
      
      for(unsigned int i = 0; i < v.edges.size(); i++)
	{
	  VID wi = v.edges[i];
	  VertInfo& w(verts[wi]);
	  if(w.index == -1) {
	    SCC_Compute(wi);
	    v.lowlink = min(v.lowlink, w.lowlink);
	  } else if(w.index > -1) {
	    v.lowlink = min(v.lowlink, w.index);
	  }
	}
      
      if(v.lowlink == v.index)
	{
	  VID wi=-1;
	  int gid = groups.size();
	  groups.push_back(vector<int>());
	  do{
	    wi = S.top(); 
	    verts[wi].group_id = gid;
	    groups.back().push_back(wi);
	    S.pop();
	  }while(wi != vi);
	}
    }
    
    // Compute the graph from the CFG.
    void SCC_Init(CFG& g)
    {
      // We're not worrying about left/right recursive.
      // We only compute the strongly connected components.
      for(int vi = 0; vi < g.nVars(); vi++)
      {
        verts.push_back(VertInfo()); 
        for(unsigned int pi = 0; pi < g.prods[vi].size(); pi++)
        {
          vector<Sym>& r(g.rules[g.prods[vi][pi].rule]);
          for(unsigned int si = 0; si < r.size(); si++)
          {
            if(r[si].isVar ())
              verts.back().edges.push_back(r[si].symID ());       
          }
        }
      }
      verts.push_back(VertInfo());
      for(int vi = 0; vi < g.nVars(); vi++)
	verts.back().edges.push_back(vi);
    }
    
  }; /* end class CFGConnect */

} /* end namespace covenant */

#endif /* __CFG_H__ */
