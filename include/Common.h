/* Common classes and types */

#ifndef __COMMON_H__
#define __COMMON_H__

#include <vector>
#include <set>
#include <avy/AvyDebug.h>

namespace covenant {

using namespace std;

  // Exception for internal errors
  class error {
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
  void remove(vector<T>& vec, size_t pos){
    assert(pos < vec.size() && "Out-of-bounds access");
    typename vector<T>::iterator it = vec.begin();
    advance(it, pos);
    vec.erase(it);
  }

  template <typename T>
  void remove_duplicates(vector<T>& vec){
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
    inline string to_string(const witness_t& witness)
    {
      stringstream ss;
      for(unsigned int ii = 0; ii < witness.size(); ii++)
      {
        ss << (char) witness[ii];
      }
      return ss.str();
    }
  } // end namespace
  
} // end namespace covenant 

#endif  /* __COMMON_H__ */
