#ifndef __CFG_PROBLEM_H__
#define __CFG_PROBLEM_H__

#include <CFG.h>

namespace covenant {

using namespace std;

class CFGProblem 
{

 public:
  class Constraint 
  {
   public:
    Constraint(const vector<int>& _vars, const CFG _lang) : 
        vars(_vars), lang(_lang){ }
    vector<int> vars;
    CFG lang;
  };

 private:
  vector<Constraint> cns;
  typedef boost::unordered_map<string, int> VarTable; 
  VarTable var_table; 
  vector<string> names;

 public:

  // Given a variable name, return the corresponding
  // string identifier.
  int getVar(string& str)
  {
    VarTable::iterator it(var_table.find(str));
    if(it != var_table.end())
    {
      return (*it).second;
    } 
    else 
    {
      int id = names.size();
      var_table[str] = id;
      names.push_back(str);
      return id;
    }
  }

  string var_name(int var)
  {
    assert(var < (int) names.size());
    return names[var];
  }

  void push(const vector<int>& vars, const CFG &lang)
  {
    cns.push_back(Constraint(vars, lang));   
  }

  Constraint& operator[](int id)
  {
    return cns[id];
  }

  int size(void) const 
  {
    return cns.size();
  }

  int nVars(void) const
  {
    return names.size();
  }

};

} /* end namespace covenant */

#endif /* __CFG_PROBLEM_H__ */
