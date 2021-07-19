#ifndef _ATOMICORBITALS_H_
#define _ATOMICORBITALS_H_

#include <vector>
using std::vector;

#include "atom.h"
/***************************************************************************************/ 
/*
 * Class define the atomic orbital for each method and get info about it
*/
/***************************************************************************************/ 
class AtomicOrbital{
 public:
  AtomicOrbital();
/***************************************************************************************/ 
   // Variables
/***************************************************************************************/ 

/***************************************************************************************/ 
   // Methods
  void setNAtom(int);
  void setElement(int);
  void setAngularMomentum(int *);
  void setIndexAO(int);

  int nAtom;
  int element;
  int angularMomentum[3];
  int indexAO;
 private:

};
/***************************************************************************************/ 
/*
 * Class define the list of all atomic orbital of the molecule 
 */
/***************************************************************************************/ 
class ListAtomicOrbitals{
 public:
  ListAtomicOrbitals();
  /***************************************************************************************/ 
  // Variables
  /***************************************************************************************/ 
  vector<AtomicOrbital> orbital;
  bool statusData;
  /***************************************************************************************/ 
  // Methods
  /***************************************************************************************/ 
  void setOrbitals(const vector<Atom> &);
 private:
  int setNumberValenceOrbitals(const int &);
  void setEachAngularMomentum(int, int* );
};



#endif // _ATOMICORBITALS_H_
