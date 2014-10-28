#ifndef __ADT_MAPSET_H__
#define __ADT_MAPSET_H__

#include <vector>
#include <set>
#include <boost/unordered_map.hpp>

namespace covenant {

using namespace std;

  // Quick class to map keys to sets
  template< typename Key, typename Value, typename Hash, typename Eq>
  class map_set {

    typedef set<Value> map_set_value_t;
    typedef boost::unordered_map< Key, map_set_value_t, Hash, Eq > map_set_t;

    typedef typename map_set_t::iterator iterator;
    typedef typename map_set_t::const_iterator const_iterator;

    map_set_t _map_set;
    
   public:
    map_set() {}
    
    ~map_set(){ _map_set.clear(); }

    map_set(const map_set&other): _map_set(other._map_set){ }
    
    void insert(Key x, Value y) {
      iterator It = _map_set.find(x);
      if (It == _map_set.end()){
        map_set_value_t s;
        s.insert(y);
        _map_set.insert(make_pair(x, s));
      }
      else{
        map_set_value_t s =  _map_set[x]; 
        s.insert(y);
        _map_set[x] =  s;
      }
    }

    set<Value> &operator[]( const Key &x) {
      return _map_set[x];
    }
        
    ostream& write(ostream& o) const {
      for(const_iterator I = _map_set.begin(), E = _map_set.end(); I!=E; ++I)
        o << (*I).first << " -> [" << (*I).second << "]" << endl;
      return o;
    }

    friend ostream& operator<<(ostream& o,  map_set m) {
      m.write(o);
      return o;
    }

  }; //end class map_set


} // end namespace covenant 

#endif  /* __ADT_MAPSET_H__ */
