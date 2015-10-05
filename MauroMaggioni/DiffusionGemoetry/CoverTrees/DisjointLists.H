#ifndef DISJOINT_LISTS_H
#define DISJOINT_LISTS_H
#include <iostream>

using namespace std;

class DisjointLists { 
private:
  int N;
  int* address;
public:
  DisjointLists(int,int*,bool = true);
  ~DisjointLists() {}
  int getNumber() const { return N; }
  int getFirst(int key) const { return address[4*key];}
  int getLast(int key) const { return address[4*key+1];}
  /* bbcrevisit
   *>> mex covertree.cpp cover.cpp DLList.cpp DisjointLists.cpp Lookup.cpp Points.cpp Vectors.cpp Sort.cpp
DisjointLists.cpp 
d:\klspectralanalysis\mariomaggioni\diffusiongemoetry\covertrees\disjointlists.h(17) : error C4716: 'DisjointLists::setFirst' : must return a value 
d:\klspectralanalysis\mariomaggioni\diffusiongemoetry\covertrees\disjointlists.h(18) : error C4716: 'DisjointLists::setLast' : must return a value 
   **/
  void setFirst(int parent,int child) { address[4*parent]=child; }
  void setLast(int parent,int child) { address[4*parent+1]=child; }
  int getPrev(int key) const { return address[4*key+2]; }
  int getNext(int key) const { return address[4*key+3]; }
  void setPrev(int key,int prev) { address[4*key+2]=prev; }
  void setNext(int key,int next) { address[4*key+3]=next; }
  void append(int parent,int key);
  void prepend(int parent,int key);
  void insertBefore(int parent,int child, int child_to_insert);
  void insertAfter(int parent,int child, int child_to_insert);
  void remove(int parent,int child);
  void clear();
  void printOn(ostream&) const;
  bool isPresent(int,int) const;
};
  
#endif
