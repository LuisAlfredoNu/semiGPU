#ifndef _ATOMICORBITALS_CPP_
#define _ATOMICORBITALS_CPP_

#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "mymath.h"

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

void AtomicOrbital::SetNAtom(int value){
  nAtom = value;
}

void AtomicOrbital::SetElement(int value){
  element = value;
}

void AtomicOrbital::SetAngularMomentum(int value[3]){
  for(int i=0;i<3;i++){
    angularMomentum[i] = value[i];
  }
}

void AtomicOrbital::SetIndexAO(int value){
  indexAO = value;
}

void AtomicOrbital::SetCoordinates(const double (&coor)[3]){
  coordinates[0] = convertAngstrom2AU(coor[0]);
  coordinates[1] = convertAngstrom2AU(coor[1]);
  coordinates[2] = convertAngstrom2AU(coor[2]);
}
/***************************************************************************************/ 
/***************************************************************************************/ 

ListAtomicOrbitals::ListAtomicOrbitals(){
  statusData = true;
}
/***************************************************************************************/ 
void ListAtomicOrbitals::SetOrbitals(const vector<Atom> &molecule){

  int moleculeSize = molecule.size();
  int totalOrbitals = 0;
  int indexAO = 0;
  
  for (int i = 0; i < moleculeSize; i++ ) {
    totalOrbitals += SetNumberValenceOrbitals(molecule[i].atomNumber);
  }

  orbital.resize(totalOrbitals);

  int angularMomentum[3];

  for (int i = 0; i < moleculeSize; i++ ) {
    if (molecule[i].atomNumber <= 2) {
      orbital[indexAO].SetNAtom(i);
      orbital[indexAO].SetElement(molecule[i].atomNumber);
      SetEachAngularMomentum(0,angularMomentum);
      orbital[indexAO].SetAngularMomentum(angularMomentum);
      orbital[indexAO].SetIndexAO(indexAO);
      orbital[indexAO].SetCoordinates(molecule[i].atomCoordinates);
      indexAO++;
      continue;
    }
    if (molecule[i].atomNumber <= 18) {
      for (int ii = 0; ii < SetNumberValenceOrbitals(molecule[i].atomNumber); ii++ ){
        orbital[indexAO].SetNAtom(i);
        orbital[indexAO].SetElement(molecule[i].atomNumber);
        SetEachAngularMomentum(ii,angularMomentum);
        orbital[indexAO].SetAngularMomentum(angularMomentum);
        orbital[indexAO].SetIndexAO(indexAO);
        orbital[indexAO].SetCoordinates(molecule[i].atomCoordinates);
        indexAO++;
      }
      continue;
    }

  }

}
/***************************************************************************************/ 
int ListAtomicOrbitals::SetNumberValenceOrbitals(const int &atomNumber){
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
void ListAtomicOrbitals::SetEachAngularMomentum(int type, int angularMomentum[3]){
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
