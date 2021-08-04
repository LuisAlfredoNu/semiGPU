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
  int nAtom;
  int element;
  int angularMomentum[3];
  int angularMomentumInt;
  int indexAO;
  double coordinates[3];
/***************************************************************************************/ 
   // Methods
/***************************************************************************************/ 
  void SetNAtom(int);
  void SetElement(int);
  void SetAngularMomentum(int *);
  void SetIndexAO(int);
  void SetCoordinates(const double (&coor)[3]);

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
  void SetOrbitals(const vector<Atom> &);
 private:
  int SetNumberValenceOrbitals(const int &);
  void SetEachAngularMomentum(int, int* );
};



#endif // _ATOMICORBITALS_H_
