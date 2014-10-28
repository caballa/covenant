#ifndef __PRODUCT_H__
#define __PRODUCT_H__

#include <boost/optional.hpp>
#include <queue>

#include <regsolver/RegSolver.h>
#include <DFA.h>

//////////////////////////////////////////////////////////////////////
// Classical textbook algorithm for intersecting regular languages
// based on the product construction.
//////////////////////////////////////////////////////////////////////

using namespace std;

namespace covenant {

template<typename EdgeSym>
class Product: public RegSolver<EdgeSym>
{    

  typedef typename RegSolver<EdgeSym>::dfa_t dfa_t;

  bool found_empty_string;

 public:

  Product(): found_empty_string(false) {}

  witness_t intersection(const vector<dfa_t> &reg_langs, 
                         unsigned int alphstart, 
                         unsigned int alphsz,
                         const int shortest_witness, 
                         bool &result)
  {      

    // Special case: all languages accept empty language
    if (!found_empty_string && allAcceptEpsilon(reg_langs))
    {
      result = true;
      return witness_t();
    }

    if (reg_langs.empty())
    {
      throw error("regular solver expects at least one regular language");
    }

    if (reg_langs.size() < 2)
    {
      result=true;
      return findShortestWitness(reg_langs[0],reg_langs[0].startState().id);
    }

    if (shortest_witness > 1)
    {
      // We generate an automata A that recognizes only strings with
      // length >= shortest_witness. Then, we intersect the input
      // automata with A.
      vector<dfa_t> reg_langs_copy;
      for(unsigned i=0; i < reg_langs.size(); i++)
      {
        reg_langs_copy.push_back(reg_langs[i]);
      }
      dfa_t fa = gen_sigma_ge_length(shortest_witness, alphstart, alphsz);
      reg_langs_copy.push_back(fa);


      boost::optional<dfa_t> tmp_inter = intersect_all(reg_langs_copy, 
                                                       alphstart, 
                                                       alphsz);
      if (tmp_inter)
      {
        result = true;
        return findShortestWitness(*tmp_inter,(*tmp_inter).startState().id);
      }
      else
      {
        // if we cannot find a witness of length >= shortest_witness
        // then we try again from scratch with the original input
        // automata before we can claim the intersection is empty.
        boost::optional<dfa_t>  inter = intersect_all(reg_langs, 
                                                      alphstart, 
                                                      alphsz);
        result = (inter ? true : false);

        if (!inter)
          return witness_t ();
        else
          return findShortestWitness(*inter,(*inter).startState().id);
      }
    }
    else
    {

        boost::optional<dfa_t>  inter = intersect_all(reg_langs, 
                                                      alphstart, 
                                                      alphsz);

        result = (inter ? true : false);

        if (!inter)
          return witness_t ();
        else
          return findShortestWitness(*inter,(*inter).startState().id);
    }
  }
    
 private:

  // PRE: sfa.size >= 2
  boost::optional<dfa_t> intersect_all(const vector<dfa_t>& reg_langs, 
                                       unsigned int alphstart, 
                                       unsigned int alphsz)
  {
    dfa_t reg_lang_0 = reg_langs[0].makeDFA (alphstart, alphsz);
    dfa_t reg_lang_1 = reg_langs[1].makeDFA (alphstart, alphsz);

    LOG("verbose", std::cout << "Automata 0 has " << reg_lang_0.nStates() << endl);
    LOG("verbose", std::cout << "Automata 1 has " << reg_lang_1.nStates() << endl);

    boost::optional<dfa_t> inter  = dfa_t::intersection(reg_lang_0, reg_lang_1);
    bool result = (inter ? true : false);
    if (result)
    {
      dfa_t prev = *inter;

      for(unsigned int i=2; i < reg_langs.size(); i++)
      {
        dfa_t reg_lang_i = reg_langs[i].makeDFA (alphstart, alphsz);

        LOG ("verbose", 
             cout << "Automata " << i << " has " << reg_lang_i.nStates() << endl);

        boost::optional<dfa_t> next = dfa_t::intersection(prev, reg_lang_i);

        result &= (next ? true : false);

        if (!result)
          break;
        
        prev = *next;
      }

      if (result)
        return boost::optional<dfa_t> (prev);
    }
    return boost::optional<dfa_t> ();
  }

  // Post: this method does not return empty witnesses even if
  // there is one.
  //
  // FIXME: it assumes that EdgeSym can be casted to int
  witness_t findShortestWitness(const dfa_t& reg_lang, unsigned int start)
  {
    // Compute shortest path
    vector<bool> visited;
    vector<int>  dist;
    vector<pair<int, int > > pi;
    for (unsigned int i=0; i< (unsigned) reg_lang.nStates(); i++)
    {
      visited.push_back(false);
      pi.push_back(make_pair(-1,-1));
      dist.push_back(-1);
    }

    queue< unsigned int> worklist;
    worklist.push(start);
    visited[start] = true;
    dist[start] = 0;
    unsigned int q=0;
    while (!worklist.empty())
    {
      q = worklist.front();
      worklist.pop();
      if (reg_lang.accepts[q] && (dist[q] > 0))
      {
        break;
      }
      for (unsigned int j=0; j < reg_lang.trans[q].size(); j++)
      {
        if (!visited[reg_lang.trans[q][j].dest])
        {
          pi[reg_lang.trans[q][j].dest] = make_pair(q,reg_lang.trans[q][j].val.symID ());
          visited[reg_lang.trans[q][j].dest] = true;
          dist[reg_lang.trans[q][j].dest] = dist[q] + 1;
          worklist.push(reg_lang.trans[q][j].dest);
        }
      } 
    } 


    // Reconstruction of the shortest path
    witness_t witness; int p;
    typename witness_t::iterator it = witness.begin();
    while ((p = pi[q].first) != -1)
    {
      it = witness.insert(it, pi[q].second);
      q = p;
    }
    if (witness.empty()) 
    {
      throw error("regular solver returned an empty witness!");
    }
    return witness; 
  }

  // Return true if for all automata its initial state is also final.
  inline bool allAcceptEpsilon(const vector<dfa_t> &reg_langs)
  {
    bool all=true;
    for (unsigned int i=0;i<reg_langs.size();i++)
    {
      all &= reg_langs[i].accepts[(reg_langs[i].startState()).id];
    }
    return all;
  }


  // Generate all the strings with length >= @length
  inline dfa_t gen_sigma_ge_length(unsigned length, 
                                   unsigned alphstart, 
                                   unsigned alphsz)
  {
    unsigned n_states = length + 1;
    dfa_t fa;
    for (unsigned i=0; i < n_states; i++)
    {
      fa.state(i);
    }
    fa.setStart(mkState(0));
    fa.accept(mkState(n_states -1 ));
    for(unsigned i=1; i < n_states; i++)
    {
      for(unsigned s=alphstart; s < alphstart+alphsz; s++)
      {
        fa.transition(mkState(i-1), 
                      EdgeSym::mkTerm(s), 
                      mkState(i));
      }
    }
    for(unsigned s=alphstart; s < alphstart+alphsz; s++)
    {
      fa.transition(mkState(n_states-1), 
                    EdgeSym::mkTerm(s), 
                    mkState(n_states-1));      
    }
    return fa;
  }

};  

} // end namespace covenant
#endif /* __PRODUCT_H__*/
