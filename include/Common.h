/* Common classes and types */

#ifndef __COMMON_H__
#define __COMMON_H__

#include <vector>
#include <set>
#include <avy/AvyDebug.h>

#include <boost/bimap/bimap.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

namespace covenant {

using namespace std;

  class TerminalFactory
  {
        
    typedef boost::bimaps::bimap<std::string, unsigned long> bimap_t;
    
    typedef bimap_t::left_iterator  left_iterator_t;
    typedef bimap_t::right_iterator right_iterator_t;

    typedef bimap_t::left_const_iterator  left_const_iterator_t;
    typedef bimap_t::right_const_iterator right_const_iterator_t;
    
    typedef boost::shared_ptr < bimap_t > bimap_ptr;

    unsigned long _next_id;
    bimap_ptr     _bimap;
    
   public:
    
    TerminalFactory(): 
        _next_id(0), _bimap( new bimap_t ())  { }
    
    TerminalFactory(unsigned long start_id): 
        _next_id(start_id), _bimap (new bimap_t ()) { }
    
    unsigned long operator[](std::string s) 
    {
      left_iterator_t it = _bimap->left.find(s);
      if (it == _bimap->left.end()) 
      {
        unsigned long res = _next_id++;
        _bimap->insert(bimap_t::value_type(s, res));
        return res;
      }
      else 
        return it->second;
    }
    
    string remap (unsigned long t) const
    {
      right_const_iterator_t it = _bimap->right.find(t);
      if (it != _bimap->right.end())
        return it->second;
      else 
        return "anonymous_" + boost::lexical_cast<string> (t);
    }
  }; 

  typedef boost::shared_ptr <TerminalFactory> TermFactory;

  // Exception for internal errors
  class error 
  {
    string _msg;
    error();
  public:
    error(string msg): _msg(msg) { }
    ostream& write(ostream& o)  {
      o << this->_msg;
      return o;
    }
    friend ostream& operator<<(ostream& o,  error e) {
      e.write(o);
      return o;
    }
  }; // class error

  // Exception for termination causes but not considered as errors
  class Exit
  {
    string _msg;
    Exit();
  public:
    Exit(string msg): _msg(msg) { }
    ostream& write(ostream& o)  {
      o << this->_msg;
      return o;
    }
    friend ostream& operator<<(ostream& o,  Exit e) {
      e.write(o);
      return o;
    }
  }; 


  template<typename T>
  inline ostream& operator<<(ostream& o, const set<T> &s) {
    o << "{";
    for(typename set<T>::iterator it= s.begin(), et= s.end(); it!=et; ++it)
      o << *it << ";" ;
    o << "}";
    return o;
  }

  template<typename T>
  inline ostream& operator<<(ostream& o, const vector<T> &v) {
    o << "[";
    for (unsigned i=0; i < v.size(); i++)
      o << v[i] << ";" ;
    o << "]";
    return o;
  }

  // Some templates to common operations on containers
  template <typename T>
  void remove(vector<T>& vec, size_t pos)
  {
    assert(pos < vec.size() && "Out-of-bounds access");
    typename vector<T>::iterator it = vec.begin();
    advance(it, pos);
    vec.erase(it);
  }

  template <typename T>
  void remove_duplicates(vector<T>& vec)
  {
    sort(vec.begin(), vec.end());
    vec.erase(unique(vec.begin(), vec.end()), vec.end());
  }

  template<typename T>
  bool intersect_unorder_ord(vector<T> unord_s, vector<T> ord_s){
    for (unsigned i=0; i < unord_s.size(); i++){
      if (binary_search(ord_s.begin(), ord_s.end(), unord_s[i]))
	return true;
    }
    return false;
  }

  template<typename T>
  bool intersect_unorder_ord(typename vector<T>::iterator  unord_begin, 
                             typename vector<T>::iterator  unord_end, 
                             vector<T>                     ord_s){
    for(typename vector<T>::iterator It = unord_begin; It!= unord_end; ++It){
      if (binary_search(ord_s.begin(), ord_s.end(), *It))
    	return true;
    }
    return false;
  }


  // To represent a witness
  // Might be worthy to have a class for it.
  typedef std::vector<int> witness_t;



  namespace witness_impl
  {
    inline string to_string(const witness_t& witness, TermFactory tfac)
    {
      if (witness.empty ()) // epsilon
        return "\"\"";

      stringstream ss;      
      for(unsigned int ii = 0; ii < witness.size(); ii++)
      {
        // ss << (char) witness[ii];
        ss << tfac->remap (witness[ii]) << " ";
      }
      return ss.str();
    }
  } // end namespace

  
} // end namespace covenant 

#endif  /* __COMMON_H__ */
