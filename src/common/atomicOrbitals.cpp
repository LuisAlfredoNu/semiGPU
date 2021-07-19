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
  angularMomentum[0]=0;
  angularMomentum[1]=0;
  angularMomentum[2]=0;
  indexAO = 0;
}

void AtomicOrbital::setNAtom(int value){
  nAtom = value;
}

void AtomicOrbital::setElement(int value){
  element = value;
}

void AtomicOrbital::setAngularMomentum(int value[3]){
  for(int i=0;i<3;i++){
    angularMomentum[i] = value[i];
  }
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

  int angularMomentum[3];

  for (int i = 0; i < moleculeSize; i++ ) {
    if (molecule[i].atomNumber <= 2) {
      orbital[indexAO].setNAtom(i);
      orbital[indexAO].setElement(molecule[i].atomNumber);
      setEachAngularMomentum(0,angularMomentum);
      orbital[indexAO].setAngularMomentum(angularMomentum);
      orbital[indexAO].setIndexAO(indexAO);
      indexAO++;
      continue;
    }
    if (molecule[i].atomNumber <= 18) {
      for (int ii = 0; ii < setNumberValenceOrbitals(molecule[i].atomNumber); ii++ ){
        orbital[indexAO].setNAtom(i);
        orbital[indexAO].setElement(molecule[i].atomNumber);
        setEachAngularMomentum(ii,angularMomentum);
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
void ListAtomicOrbitals::setEachAngularMomentum(int type, int angularMomentum[3]){
    if (type == 0) {
      angularMomentum[0] = 0;
      angularMomentum[1] = 0;
      angularMomentum[2] = 0;
      return; 
    } else if (type == 1) {
      angularMomentum[0] = 0;
      angularMomentum[1] = 0;
      angularMomentum[2] = 1;
      return; 
    } else if (type == 2) {
      angularMomentum[0] = 0;
      angularMomentum[1] = 1;
      angularMomentum[2] = 0;
      return; 
    } else if (type == 3) {
      angularMomentum[0] = 1;
      angularMomentum[1] = 0;
      angularMomentum[2] = 0;
      return; 
    }
}
#endif // _ATOMICORBITALS_CPP_
