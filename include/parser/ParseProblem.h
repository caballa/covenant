#ifndef __PARSE_CFG_PROBLEM_H__
#define __PARSE_CFG_PROBLEM_H__

#include <CFGProblem.h>
#include <parser/ParseCFG.h>

// Parses a specification of language constraints.
// Hampi format is: (elem (##VARS##) (##CFL##))
// Otherwise the format is (##CFL##)

using namespace std;

namespace covenant {

template<class In>
void consume(In& input, const char* str)
{
  while(*str != '\0')
  {
    input.chomp(*str);
    str++;
  }
}

// A variable identifier.
// For now, we assume all identifiers are entirely alpha.
template<class In>
string parse_ident(In& input)
{
  string name;
  consume_blanks(input);

  while(!input.empty() && isalpha(input.peek()))
    name.push_back(input.pop());

  if (!( name.size() > 0))
    throw error("make sure CFG is enclosed between parenthesis");

  return name;
}

template<class In>
void parse_constraint(CFGProblem &p, In& input)
{
  vector<int> vars;
  consume_blanks(input);

  if (input.empty())
    throw error("make sure CFG is enclosed between parenthesis");
  
#ifdef HAMPI_FORMAT
  input.chomp('('); // open constraint
  consume_blanks(input);
  consume(input, "elem");

  consume_blanks(input);
  input.chomp('('); // open list
  consume_blanks(input);
  char c = input.peek();
  while(c != ')')
  {
    // Read another var name
    string varname = parse_ident(input);
    int varid = p.getVar(varname);
    vars.push_back(varid);

    // Progress input
    consume_blanks(input);
    if (input.empty())
      throw error("make sure CFG is enclosed between parenthesis");

    c = input.peek();
  }
  input.chomp(); // close list
  consume_blanks(input);
#endif 

  input.chomp('('); // open cfg
  consume_blanks(input);
  CFG lang = parse_cfg(input);
  input.chomp(')'); // close cfg
  consume_blanks(input);

#ifdef HAMPI_FORMAT
  input.chomp(')'); // close constraint
#else
  string varname("x");
  int varid = p.getVar(varname);
  vars.push_back(varid);
#endif 
  p.push(vars, lang);
}

template<class In>
void parse_problem(CFGProblem& p, In& input)
{
  consume_blanks(input);
  while(!input.empty())
  {
    char c = input.peek();
    // gross fix if the last constraint does not finish with some blank
    if (input.peek() == ')') break;

    if(c == ';')
    {
      // Comment line.
      consume_comment(input);
    } 
    else 
    {
      if (input.peek() != '(')
        throw error("make sure CFG is enclosed between parenthesis");

      consume_blanks(input);
      parse_constraint(p, input);
      consume_blanks(input);
    }
  }
}

} // end namespace covenant 

#endif /*__PARSE_CFG_PROBLEM_H__*/
