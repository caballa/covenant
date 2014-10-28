#ifndef __CFG_SOLVER_H__
#define __CFG_SOLVER_H__

#include <fstream>

#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include <Common.h>
#include <CFG.h>
#include <DFA.h>
#include <CFGProblem.h>

#include <abstraction/Abstract.h>
#include <abstraction/SigmaStarAbstract.h>
#include <abstraction/CycleBreakingAbstract.h>

#include <refinement/EpsRefine.h>
#include <refinement/GreedyEpsRefine.h>
#include <refinement/MaxEpsRefine.h>

#include <regsolver/RegSolver.h>
#include <regsolver/Product.h>

namespace covenant{

using namespace std;

enum STATUS { SAT, UNSAT, UNKNOWN};

namespace witness_impl
{
  template<typename EdgeSym>
  inline DFA<EdgeSym> to_automata(const witness_t &witness)
  {
    if (witness.empty())
      throw error("regular solver returned an empty witness!");

    DFA<EdgeSym> dfa;
    State start, final;
    State prev , curr;
    unsigned n = witness.size();
    for(unsigned i=0; i<= witness.size(); i++)
    {
      curr = dfa.state(i);
      if (i == 0)
        start = curr;
      else
      {
        if (i == n)
        {
          final = curr;
        }
        dfa.transition(prev, EdgeSym::mkTerm(witness[i-1]), curr);
      }
      prev = curr;
    }
    dfa.setStart(start);
    dfa.accept(final);
    return dfa;
  }
} // end namespace

template<typename EdgeSym>
class Solver: public boost::noncopyable {

 public:

  typedef RegSolver<EdgeSym>             reg_solver_t;
  typedef DFA<EdgeSym>                   dfa_t;
  typedef CFGDigest<EdgeSym>             cfg_digest_t;
  typedef CondEpsGen<EdgeSym>            cond_eps_gen_t;
  typedef SigmaStarAbstract<EdgeSym>     sigma_start_abs_t;
  typedef CycleBreakingAbstract<EdgeSym> cycle_breaking_abs_t;
  typedef EpsRefine<EdgeSym>             eps_refine_t;
  typedef MaxEpsRefine<EdgeSym>          max_eps_refine_t;
  typedef GreedyEpsRefine<EdgeSym>       greedy_eps_refine_t;

  struct Options
  {
    Options(int o1, GeneralizationMethod o2, AbstractMethod o3,
            bool o4=false, int o5=1, int o6=-1, unsigned o7=1):
        max_cegar_iter(o1),
        gen(o2), 
        abs(o3),
        is_dot_enabled(o4),
        shortest_witness(o5),
        freq_incr_witness(o6),
        incr_witness(o7){ } 
    int max_cegar_iter;
    GeneralizationMethod gen;
    AbstractMethod abs;
    bool is_dot_enabled;
    int  shortest_witness;
    int  freq_incr_witness;
    unsigned incr_witness;
  };


 private:

  // Contains the vector of CFGs to be intersected
  CFGProblem problem; 

  // Regular solver 
  reg_solver_t *solver;

  // Language digests
  vector<CFG> cfgs;
  vector<cfg_digest_t> digests;
  
  // Regular approximations
  vector<dfa_t> reg_langs;
  vector<bool> exact;

  // Solver options
  Options opts;

  // Logging/debugging purposes
  ofstream refine_log; 
  ofstream abstract_log; 
  ofstream proof_log; 

  ////
  // Methods for implementing the CEGAR loop
  ////

  // Perform initial abstraction: convert each cfg to a finite
  // automata that overapproximates the language of the context-free
  // grammar.  If is_regular[i] is true then the conversion does not
  // lose precision.
  vector<dfa_t> abstraction(const vector<CFG>  &cfgs, 
                            const vector<bool> &is_regular, 
                            AbstractMethod abs)
  { 
    vector<dfa_t> fas;
    fas.reserve(cfgs.size());
    for(unsigned int i = 0 ; i < cfgs.size() ; i++)
    {
      dfa_t fa;
      switch (abs)
      {
        case SIGMA_STAR: 
          {
            sigma_start_abs_t alpha; 
            fa = alpha.do_abstraction(cfgs[i], is_regular[i]);
          }
          break;
        default: 
          {
            cycle_breaking_abs_t alpha;
            fa = alpha.do_abstraction(cfgs[i], is_regular[i]);
          }
      }

      if (opts.is_dot_enabled)
      {
        fa.print_dot(abstract_log, 
                     std::string ("Initial Regular Approximation of CFG ") + 
                     boost::lexical_cast<string> (i));          
      }

      LOG ("solver", cout << "Approximated Finite Automata:\n" << fa);
      fas.push_back(fa);
    }
    return fas;
  }
  
  // Produce the next regular approximations for the cegar loop.
  bool refine(witness_t& witness, 
              GeneralizationMethod gen, 
              bool REFINE_ONLY_FIRST) 
  {

    LOG ("solver", 
         cout << "Witness found: " << witness_impl::to_string(witness) << "\n");
    
    if (opts.is_dot_enabled)
    {
      dfa_t cex = witness_impl::to_automata<EdgeSym>(witness);
      cex.print_dot(refine_log, "Cex ") ;
    }

    
    // Choose here the generalization method.
    eps_refine_t* r = NULL;
    if(gen == MAX_GEN)
      r = new max_eps_refine_t(cfgs[0].alphstart, cfgs[0].alphsz);
    else 
      r = new greedy_eps_refine_t();
    
    bool refined = false;
    for(unsigned int li = 0; li < digests.size(); li++)
    {
      // L_i is regular, so w \in L_i
      if(exact[li]) 
      {
        LOG ("verbose" , std::cout << "Automata " << li << " is exact\n");
        continue;
      }
      
      // First, we check if the witness is in the language.
      cond_eps_gen_t gen(digests[li], witness);
      if(!gen.is_empty())
      {
        LOG ("verbose", std::cout << "Cex is accepted by CFG " << li << "\n");
        continue;
      }
      
      // Otherwise, we refine it.
      dfa_t gen_lang(r->refine(gen));
      
// #ifndef REFINE_DIFF
//       // Haven't currently extended solvers to handle complemented
//       // languages.
//       gen_lang.eps_elim();
//       dfa_t ref_complement(complement(makeDFA(gen_lang, 
//                                                       cfgs[li].alphstart, 
//                                                       cfgs[li].alphsz)));
      
//       // We can either just add ref_complement to the set of regular
//       // languages, or intersect it with an approximation.
//       reg_langs.push_back(ref_complement);
// #else
      gen_lang.eps_elim();

      if (opts.is_dot_enabled)
      {
        gen_lang.print_dot(refine_log, 
                           string("Generalized cex wrt CFG ") + 
                           boost::lexical_cast<string> (li)) ;
      }

      boost::optional<dfa_t> refined_fa = dfa_t::difference (reg_langs[li], 
                                                             gen_lang, 
                                                             cfgs[li].alphstart, 
                                                             cfgs[li].alphsz);
      
      if (!refined_fa)
      {
        throw error ("refinement did not make any progress");
      }
      
      reg_langs[li] = *refined_fa;
 
      LOG ("solver", 
           cout << "New refined regular approximation:\n" << *refined_fa << "\n");
           

//#endif
      refined = true;
      LOG ("verbose", std::cout << "Automata " << li << " refined.\n");
      if(REFINE_ONLY_FIRST)
        break;
    }
    
    delete r;
    
    return refined;
  }

 public:

  // Constructor of the class
  Solver(CFGProblem problem_, reg_solver_t *s, const Options opts_): 
      problem(problem_), solver(s), opts(opts_)
  {  
    if (opts.is_dot_enabled)
    {
      refine_log.open("refinements.dot");
      abstract_log.open("abstractions.dot");
      proof_log.open("proof.dot");
    }
  }
    
  // Destructor of the class
  ~Solver()
  { 
    if (opts.is_dot_enabled)
    {
      refine_log.close();
      abstract_log.close();
      proof_log.close();
    }
  }

  // The solver implements a CEGAR loop. Since the intersection of CFGs
  // is an undecidable problem the solver may return "yes"
  // (ie. intersection not empty), "no" (i.e, intersection is empty) or
  // run forever.
  STATUS solve()
  {
    
    witness_t witness;
    const bool REFINE_ONLY_FIRST = false;

    int iter=1;
    cfgs.clear();
    
    exact.reserve(problem.size());
    bool all_cfgs_are_reg=true;
    for(int ii = 0; ii < problem.size(); ii++)
    {
      CFGProblem::Constraint cn(problem[ii]);
      // This must be done before we normalize the cfg.
      const bool is_cfg_reg = cn.lang.isRegularGrammar();
      exact.push_back(is_cfg_reg);
      all_cfgs_are_reg &= is_cfg_reg;
      cn.lang.normalize();
      LOG("solver", 
          cout << "After normalization:" << endl <<  cn.lang << endl;
          cout << "===================================================\n");
      cfgs.push_back(cn.lang);
    }
    
    unsigned int alphstart;
    unsigned int alphsz;
    shrink_alphabet(cfgs, alphstart, alphsz);
#ifdef SANITY_CHECKS
    // Ensure all cfgs have the same alphabet
    for (unsigned int i=0;i<cfgs.size();i++){
      const bool same_alph =  ( cfgs[i].alphstart == (int) alphstart &&
                                cfgs[i].alphsz    == (int) alphsz );
      if (!same_alph)
        throw error("solver expects all the CFGs with same alphabet");
    }
#endif

    reg_langs = abstraction(cfgs, exact, opts.abs);
    
    LOG ("solver", 
         cout << "Initial regular approximations done." << endl;
         cout << "===================================================\n");
    
    for(unsigned int li = 0; li < cfgs.size(); li++)
      digests.push_back(cfg_digest_t(cfgs[li]));

    LOG ("solver", 
         cout << "CFG digests computed." << endl;
         cout << "===================================================\n");

    ////
    // CEGAR loop
    ////
    int freq_incr_witness= opts.freq_incr_witness;
    while (true){
      if (opts.max_cegar_iter >= 0 && iter == opts.max_cegar_iter) 
        return UNKNOWN;

      LOG( "verbose" , 
           cout << "Checking intersection of regular approximations ...\n");

      // solver between regular languages
      // pre: automata do no have epsilon transitions
      bool sat_query;

      // This increases the length of the witnesses so we don't get
      // stuck

      int shortest_witness = opts.shortest_witness; 

      if (freq_incr_witness == 0)
      {
        shortest_witness += opts.incr_witness;
        LOG( "verbose", cout << "Searching for witnesses of length " 
                             << shortest_witness << endl); 
        freq_incr_witness = opts.freq_incr_witness;
      }
      else if (freq_incr_witness > 0)
      {
        freq_incr_witness--;
      }

      witness = solver->intersection(reg_langs, 
                                     alphstart, 
                                     alphsz,
                                     shortest_witness,
                                     sat_query);

      LOG("solver", 
          cout << "Intersection done.\n";
          cout << "===================================================\n");

      if (!sat_query)
      {
        LOG( "verbose", cout << "Intersection of regular approximations is empty!\n");
        cout << "Finished after " << iter << " cegar iterations.\n";
        if (opts.is_dot_enabled)
        {
          for(unsigned ii=0;ii<cfgs.size();ii++)
          {
            reg_langs[ii].print_dot(proof_log, 
                                    string ("Final Regular Approximation of CFG ") + 
                                    boost::lexical_cast<string>(ii));
          }
        }
        return UNSAT;
      }
      else if (all_cfgs_are_reg)
      {
        cout << witness_impl::to_string(witness) << endl;
        cout << "Found a real cex after " << iter << " iterations.\n";
        return SAT;
      }
      else
      {
        if(!refine(witness, opts.gen, REFINE_ONLY_FIRST))
        {
          cout << witness_impl::to_string(witness) << endl;
          cout << "Found a real cex after " << iter << " iterations.\n";
          return SAT;
        }
      }

      LOG ("solver", 
           cout << "Refinement done.\n";
           cout << "===================================================\n");

      iter++;
    }
  }

 private:
    
  // Count for all terminal symbols to make smaller the alphabet of the
  // grammars. This can have a significant impact later on for some
  // finite automata operations.
  // POST: modify the attributes alphstart and alphsz from the
  // grammars but the production rules are not modified.
  inline void shrink_alphabet(vector<CFG> &cfgs, 
                              unsigned &alphstart, 
                              unsigned &alphsz)
  {
    if (cfgs.empty())
      throw error("solver expects at least one CFG");
    
    int min, max;
    cfgs[0].get_alphabet_min_and_max(min, max);    
    for (unsigned int i=1; i<cfgs.size(); i++)
    {
      int g_min, g_max;
      cfgs[i].get_alphabet_min_and_max(g_min, g_max);
      min = std::min(min, g_min);
      max = std::max(max, g_max);
    }
    alphstart = min;
    alphsz    = (max - min) + 1;
    // Shrink (if possible) all grammar's alphabets
    for (unsigned int i=0; i<cfgs.size(); i++)
    {
      cfgs[i].alphstart = alphstart;
      cfgs[i].alphsz    = alphsz;
    }
  }

};

} // end namespace covenant
#endif /*__CFG_SOLVER_H__*/

