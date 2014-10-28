/* 
 * This file is from the revenant tool (http://ww2.cs.mu.oz.au/~ggange/revenant/)
 */

#ifndef __ADT_SPARSE_SET_H__
#define __ADT_SPARSE_SET_H__

// Data structure for maintaining monotonically increasing
// sets of elements within a bounded range.
// Permits constant time lookup, insertion and iteration. 

// U indicates whether to use FAST_FSET. (Slightly more expensive
// insertion, balanced by slightly cheaper lookups.)
// In most cases, it's small enough that it won't make any difference.
#include <cassert>
#include <cstdlib>

template<int U = 0>
class SparseSet {

 protected:
  unsigned int  dom;
  unsigned int* sparse;
  unsigned int* dense;
  unsigned int  members;

 public:

  SparseSet(void) : dom(0), sparse(NULL), dense(NULL), members( 0 )
  { }
  
  SparseSet(unsigned int size) : 
    dom(size), sparse(new unsigned int[size]), 
    dense(new unsigned int[size]), members( 0 )
  {
    if( U&1 )
    {
      assert( members == 0 );
      for( unsigned int i = 0; i < dom; i++ )
      {
        sparse[i] = i;
        dense[i] = i;
      }
    }
  }
  
  ~SparseSet() 
  {
    if( sparse )
      delete[] sparse;
    if( dense )
      delete[] dense;  
  }
  
  bool elem(unsigned int value) const 
  {
    if( U&1 )
    {
      return (sparse[value] < members);
    } 
    else 
    {
      unsigned int a = sparse[value];
      if( a < members && dense[a] == value )
        return true;
      return false;
    }
  }
   
  bool elemLim(unsigned int lim, unsigned int el)
  {
    return (sparse[el] < lim) && elem(el);
  }
  
  virtual bool insert(unsigned int value)
  {
    if( U&1 )
    {
      unsigned int old_dense = sparse[value];
      unsigned int lost_val = dense[members];
      
      sparse[value] = members;
      dense[members] = value;
      
      sparse[lost_val] = old_dense;
      dense[old_dense] = lost_val;
    } 
    else 
    {
      assert( !elem(value) );
      
      sparse[value] = members;
      dense[members] = value;
    }
    members++;
    return true;
  }
  
   void clear(void) 
  {
     members = 0;
   }

   unsigned int pos(unsigned int val) const
   {
     assert( elem(val) );
     return sparse[val];
   }
  
   unsigned int operator [] (unsigned int index) 
  {
     assert(index < members);
     return dense[index];
   }
    
   void growTo(unsigned int sz)
   {
     if( sz > dom )
     {
       sparse = (unsigned int*) realloc(sparse,sizeof(unsigned int)*sz); 
       dense = (unsigned int*) realloc(dense,sizeof(unsigned int)*sz);
       
       if( U&1 )
       {
         assert( members == 0 );
         for(; dom < sz; dom++ )
         {
           sparse[dom] = dom;
           dense[dom] = dom;
         }
       }
     }
   }
  
  unsigned int size(void) 
  {
    return members;
  }
   
  unsigned int domain(void) 
  {
      return dom;
  }
       
};

#endif
