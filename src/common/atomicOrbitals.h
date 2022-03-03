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
  int indexAtom;
  int element;
  int angularMomentum[3];
  int angularMomentumInt;
  int indexAO;
  double coordinates[3];
/***************************************************************************************/ 
   // Methods
/***************************************************************************************/ 
  void SetIndexAtom(int);
  void SetElement(int);
  void SetAngularMomentum(int *);
  void SetIndexAO(int);
  void SetCoordinates(const double (&coor)[3]);
  #pragma acc routine seq
  int GetCoreCharge() const;
  #pragma acc routine seq
  int GetAOsSize() const;

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
  ~ListAtomicOrbitals();
  /***************************************************************************************/ 
  // Variables
  /***************************************************************************************/ 
  AtomicOrbital* orbital;
  bool statusData;
  /***************************************************************************************/ 
  // Methods
  /***************************************************************************************/ 
  void SetOrbitals(const vector<Atom> &);
  int GetOrbital4NextAtom(const int&,const vector<AtomicOrbital>&);
  #pragma acc routine seq
  size_t size() const;
 private:
  /***************************************************************************************/ 
  // Variables
  /***************************************************************************************/ 
  size_t totalOrbitals_;
  /***************************************************************************************/ 
  // Methods
  /***************************************************************************************/ 
  int SetNumberValenceOrbitals(const int &);
  void SetEachAngularMomentum(int, int* );
  /***************************************************************************************/ 
  // OpenAcc
  /***************************************************************************************/ 
  void To_device();
  void From_device();
  void Update_host();
  void Update_device();

};

#endif // _ATOMICORBITALS_H_
