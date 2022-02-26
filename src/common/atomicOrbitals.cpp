#ifndef _ATOMICORBITALS_CPP_
#define _ATOMICORBITALS_CPP_

#include "mymath.h"

#include "atomicOrbitals.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
AtomicOrbital::AtomicOrbital(){
  indexAtom = 0;
  element = 0;
  angularMomentum[0]=0;
  angularMomentum[1]=0;
  angularMomentum[2]=0;
  indexAO = 0;
}
/***************************************************************************************/ 
void AtomicOrbital::SetIndexAtom(int value){
  indexAtom = value;
}
/***************************************************************************************/ 
void AtomicOrbital::SetElement(int value){
  element = value;
}
/***************************************************************************************/ 
void AtomicOrbital::SetAngularMomentum(int value[3]){
  angularMomentumInt = 0;
  for(int i=0;i<3;i++){
    angularMomentum[i] = value[i];
    if (value[i] == 1) {
      angularMomentumInt = i + 1;
    }
  }
}
/***************************************************************************************/ 
void AtomicOrbital::SetIndexAO(int value){
  indexAO = value;
}
/***************************************************************************************/ 
void AtomicOrbital::SetCoordinates(const double (&coor)[3]){
  coordinates[0] = convertAngstrom2AU(coor[0]);
  coordinates[1] = convertAngstrom2AU(coor[1]);
  coordinates[2] = convertAngstrom2AU(coor[2]);
}
/***************************************************************************************/ 
int AtomicOrbital::GetCoreCharge() const {
  if (element < 2) {
    return element;
  }else if (element < 11 ) {
    return element - 2 ;
  }else if (element < 19) {
    return element - 10;
  }
  return -1;
}
/***************************************************************************************/ 
int AtomicOrbital::GetAOsSize() const {
  if (element < 2) {
    return 1;
  }else if (element < 11 ) {
    return 4 ;
  }
  return -1;
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

  totalOrbitals_ = totalOrbitals;

  orbital = new AtomicOrbital[totalOrbitals];

  int angularMomentum[3];

  for (int i = 0; i < moleculeSize; i++ ) {
    if (molecule[i].atomNumber <= 2) {
      orbital[indexAO].SetIndexAtom(i);
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
        orbital[indexAO].SetIndexAtom(i);
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
      angularMomentum[0] = 1;
      angularMomentum[1] = 0;
      angularMomentum[2] = 0;
      return; 
    } else if (type == 2) {
      angularMomentum[0] = 0;
      angularMomentum[1] = 1;
      angularMomentum[2] = 0;
      return; 
    } else if (type == 3) {
      angularMomentum[0] = 0;
      angularMomentum[1] = 0;
      angularMomentum[2] = 1;
      return; 
    }
}
/***************************************************************************************/ 
int ListAtomicOrbitals::GetOrbital4NextAtom(const int& actualOrbital,const vector<AtomicOrbital>& AOs){
  int actualAtom = AOs[actualOrbital].indexAtom;
  return  actualOrbital + SetNumberValenceOrbitals(actualAtom);
}
/***************************************************************************************/ 
size_t ListAtomicOrbitals::size() const {
  return totalOrbitals_;
}
#endif // _ATOMICORBITALS_CPP_
