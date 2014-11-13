#ifndef __PARSE_CFG_H__
#define __PARSE_CFG_H__

#include <cctype>

#include <Common.h>
#include <CFG.h>

// Simple parser for specifying  CFG.

// Interfaces
//============
// In:
// char peek(void)    -- returns the first character of input.
// void chomp(void)   -- consumes the first character of input.
// void chomp(char c) -- consumes the first character of input,
//                       fails if first character is not c.
// char pop(void)     -- consumes and returns the first character of input.
// bool empty(void)   -- checks if the input stream is empty.
//                                 

using namespace std;

namespace covenant {

typedef boost::unordered_map<string, Sym> NonTerminalMap;

template<class In>
void consume_blanks(In& input)
{
  while (!input.empty() && 
         (input.peek() == ' ' || input.peek() == '\n' || input.peek() == '\t'))
  {
    input.chomp();
  }
}

template<class In>
void consume_comment(In& input)
{
  while(!input.empty() && input.peek() != '\n')
    input.chomp();
  input.chomp();
  consume_blanks(input);
}

template<class In, class TerminalFactory>
CFG parse_cfg(In& input, TerminalFactory tfac)
{
  CFG g(tfac);

  NonTerminalMap mapping;
  parse_rules(input, g, mapping);
  if(input.empty()) return g;
  consume_blanks(input);
  // the special case if input.peek() == ')' is in case we read the
  // format (elem (##VARS##) (##CFL##))
  if (input.empty() || input.peek() == ')') 
    return g;
  
  throw error("during the parsing of the cfg");
}

template<class In>
void parse_rules(In& input, CFG &g, NonTerminalMap &mapping)
{

  char c = input.peek();
  parse_rule(input, g, mapping);
  
  if(input.empty())
    return;
  
  // Check if we continue parsing
  c = input.peek();
  switch(c)
  {
    case ';':
      // Another rule
      input.chomp();
      consume_blanks(input);
      parse_rules(input, g, mapping);
    default:
      break;
  }
  return;
}

// End of symbol
bool EOfSymbol(char c)
{
  return (c == ' ' || c == '-' || c == ']' || c == ',');
}

template<class In>
Sym parse_symbol(In& input, CFG &g, NonTerminalMap &mapping)
{
  char c;
  string name;
  while(!input.empty() && !EOfSymbol(c = input.peek()))
  {
    input.chomp();
    name += c;
  }

  LOG ("parser" , cout << "Read symbol: " << name << endl);

  NonTerminalMap::iterator it = mapping.find(name);
  if (it != mapping.end())
    return (*it).second;

  Sym s( g.newVar() );
  mapping.insert(make_pair(name,s));
  return s;
}

template<class In>
int parse_hex(In& input)
{ 
  // The format is: "\xnn"  where nn are exactly two hexadecimal digits.
  
  int val = 0;
  char c;
  while(!input.empty() && isxdigit(c = input.peek()))
  {
    input.chomp();
    if (isdigit(c))
      val = 16*val + (c - '0');
    else
      val = 16*val + (tolower(c) -'a' + 10);
  }
  return val;
}


template<class In>
std::string parse_string(In& input)
{ 
  std::string str;
  char c;
  while(!input.empty() && ( (c = input.peek()) != '\"'))
  {
    input.chomp ();
    str.push_back (c);
  }
  if (!input.empty ())
    input.chomp ();

  return str;
}

template<class In>
void parse_rule(In& input, CFG &g, NonTerminalMap &mapping)
{
  LOG ("parser" , 
       cout << "String to consume: " << input.str.substr(input.idx) << endl) ;

  Sym s = parse_symbol(input, g, mapping);
  consume_blanks(input);
  input.chomp('-');
  input.chomp('>');
  consume_blanks(input);
  input.chomp('[');
  parse_rule_rhs(input,g,s,mapping);
}

// End of one production
bool EOfProduction(char c)
{
  return c == ']';
}

template<class In>
void parse_rule_rhs(In& input, CFG &g, Sym lhs, NonTerminalMap &mapping)
{
  vector<Sym> ss;
  char c;

  c = input.peek();
  if (c == ']')
  {
    input.chomp();
    consume_blanks(input);
    // epsilon transition
    g.prod(lhs, Rule::E (g.getTermFactory ()));
    return;
  }

  while(!input.empty() && (!EOfProduction(c = input.peek())))
  {
    LOG ("parser", 
         cout << "String to consume: " << input.str.substr(input.idx)
         << endl);

    switch(c)
    {
      case ' ':
        input.chomp();
        break;
      case ',':
        {
	input.chomp();
	Rule R (g.getTermFactory ());
	for(unsigned i=0; i<ss.size(); i++)
          {
	  R << ss[i];
	}
	g.prod(lhs, R);
	ss.clear();
        }
        break;
      case '\"':
        {
	input.chomp();
          std::string str = parse_string(input);
	Sym term = Sym::mkTerm( g.getTermFactory()->operator[](str));
	ss.push_back(term);
        }
        break;
      case '\\':
        {
	input.chomp();
	input.chomp('x'); 
	int val = parse_hex(input);
	Sym term = Sym::mkTerm(val);
	ss.push_back(term);
        }
        break;
      default:
        ss.push_back(parse_symbol(input, g,mapping));
    } // end switch
  } // end while  
  
  LOG ("parser" , 
       cout << "String to consume: " << input.str.substr(input.idx)
       << endl) ;

  // last rule
  input.chomp(']');
  consume_blanks(input);
  Rule R (g.getTermFactory ());
  for(unsigned i=0; i<ss.size(); i++)
  {
    R << ss[i];
  }
  g.prod(lhs, R);
  ss.clear();
}

} // end namespace covenant 

#endif /* __PARSE_CFG_H__  */

