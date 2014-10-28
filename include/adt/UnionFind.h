#ifndef __ADT_UNION_FIND_H__ 
#define __ADT_UNION_FIND_H__ 

#include <vector>

namespace covenant {

using namespace std;

  // Quick class for an union-find
  class UnionFind {
   public:
    UnionFind(int N){
      // Each class maps to itself.
      for(int ii = 0; ii < N; ii++)
        mapping.push_back(ii);
    }
    
    int lookup(int id){
      if(mapping[id] == id)
        return id;
      int ret = lookup(mapping[id]);
      mapping[id] = ret;
      return ret;
    }

    int set(int src, int dest){
      int edest = lookup(dest);
      int esrc = lookup(src);
      mapping[esrc] = edest;
      return edest;
    }

    unsigned size(){
      return mapping.size();
    }

   private:
    vector<int> mapping;

  }; /* end class UnionFind*/

} // end namespace covenant 

#endif  /* __ADT_UNION_FIND_H__ */
