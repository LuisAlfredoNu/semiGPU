#ifndef _ATOMICORBITALS_CPP_
#define _ATOMICORBITALS_CPP_

#include <iostream>
using std::cout;
#include <vector>
using std::vector;

#include "atomicOrbitals.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

AtomicOrbital::AtomicOrbital(){
  nAtom = 0;
  element = 0;
  angularMomentum = {0,0,0};
  indexAO = 0;
}

void AtomicOrbital::setNAtom(int value){
  nAtom = value;
}

void AtomicOrbital::setElement(int value){
  element = value;
}

void AtomicOrbital::setAngularMomentum(vector<int> value){
  angularMomentum = value;
}

void AtomicOrbital::setIndexAO(int value){
  indexAO = value;
}
/***************************************************************************************/ 
/***************************************************************************************/ 

ListAtomicOrbitals::ListAtomicOrbitals(){
  statusData = true;
}
/***************************************************************************************/ 
void ListAtomicOrbitals::setOrbitals(const vector<Atom> &molecule){

  int moleculeSize = molecule.size();
  int totalOrbitals = 0;
  int indexAO = 0;
  
  for (int i = 0; i < moleculeSize; i++ ) {
    totalOrbitals += setNumberValenceOrbitals(molecule[i].atomNumber);
  }

  orbital.resize(totalOrbitals);

  for (int i = 0; i < moleculeSize; i++ ) {
    if (molecule[i].atomNumber <= 2) {
      orbital[indexAO].setNAtom(i);
      orbital[indexAO].setElement(molecule[i].atomNumber);
      int angularMomentum[3] = setEachAngularMomentum(0);
      orbital[indexAO].setAngularMomentum(angularMomentum);
      orbital[indexAO].setIndexAO(indexAO);
      indexAO++;
      continue;
    }
    if (molecule[i].atomNumber <= 18) {
      for (int ii = 0; ii < setNumberValenceOrbitals(molecule[i].atomNumber); ii++ ){
        orbital[indexAO].setNAtom(i);
        orbital[indexAO].setElement(molecule[i].atomNumber);
        vector<int> angularMomentum = setEachAngularMomentum(ii);
        orbital[indexAO].setAngularMomentum(angularMomentum);
        orbital[indexAO].setIndexAO(indexAO);
        indexAO++;
      }
      continue;
    }

  }

  cout << std::endl;
  cout << "Memory for orbital " << sizeof(AtomicOrbital) * orbital.size() << std::endl;
  cout << "Size memory with molecule of 1000 atoms " << sizeof(AtomicOrbital) * 32 << std::endl;


}
/***************************************************************************************/ 
int ListAtomicOrbitals::setNumberValenceOrbitals(const int &atomNumber){
  if (atomNumber <= 2) {
    return 1;
  } else if (atomNumber <= 18) {
    return 4;
  } else if (atomNumber > 18) {
    statusData = false;
  }
  return -1;
}
/***************************************************************************************/ 
vector<int> ListAtomicOrbitals::setEachAngularMomentum(int type){
    vector<int> angularMomentum = {0,0,0};
    if (type == 0) {
      return angularMomentum;
    } else if (type == 1) {
      angularMomentum = {0,0,1};
      return angularMomentum;
    } else if (type == 2) {
      angularMomentum = {0,1,0};
      return angularMomentum;
    } else if (type == 3) {
      angularMomentum = {1,0,0};
      return angularMomentum;
    }
}
#endif // _ATOMICORBITALS_CPP_
