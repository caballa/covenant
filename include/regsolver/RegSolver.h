#ifndef __REGULAR_SOLVER_H__
#define __REGULAR_SOLVER_H__

#include <Common.h>
#include <CFG.h>
#include <DFA.h>

namespace covenant {

using namespace std;

// Here the different regular solvers
enum RegSolverName { PRODUCT };

template<typename EdgeSym>
class RegSolver 
{    

 public:

  typedef DFA<EdgeSym> dfa_t;

 public:

  virtual ~RegSolver () { }
  virtual witness_t intersection(const vector<dfa_t>& reg_langs, 
                                 unsigned int alphstart, 
                                 unsigned int alphsz,
                                 const int shortest_witness, 
                                 bool &result) = 0;

};   // end RegSolver

} // end namespace covenant

#endif /* __REGULAR_SOLVER_H__*/
