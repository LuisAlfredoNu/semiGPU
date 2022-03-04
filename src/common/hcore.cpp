#ifndef _HCORE_CPP_
#define _HCORE_CPP_

#include <iostream>
using std::cout;
using std::endl;



#include "atomicOrbitals.h"
#include "twocenterintegral.h"
#include "mymemory.h"

#include "hcore.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
Hcore::Hcore(const MNDOparameter& MNDOpara,const ListAtomicOrbitals& infoAOs,\
             double**** all2CenInt) : BaseMatrix(infoAOs.size()){
  parameter_ = &MNDOpara;
  infoAOs_ = &infoAOs;
  all2CenterIntegral_ = all2CenInt;

  overlap_ = new Overlap();
  To_device();
}
/***************************************************************************************/ 
double Hcore::ComputeElementMatrix(const size_t &i,const size_t &j){
  // Diagonal elements
  if ( i == j ) {
    return ComputeDiagonal(infoAOs_->orbital[i]);
  }
  // Non diagonal element, same atom
  if (infoAOs_->orbital[i].indexAtom == infoAOs_->orbital[j].indexAtom) {
    return ComputeNonDiagonalSameAtom(infoAOs_->orbital[i],infoAOs_->orbital[j]);
  }else{
    return ComputeNonDiagonalDiffAtom(infoAOs_->orbital[i],infoAOs_->orbital[j]);
  }
}
/***************************************************************************************/ 
double Hcore::ComputeNonDiagonalDiffAtom(const AtomicOrbital& orbitalA,\
       const AtomicOrbital& orbitalB){

  double overlapValue = overlap_->ComputeOverlap(orbitalA,orbitalB);
  //double overlapValue = overlap_->ComputeDummy(orbitalA,orbitalB);
  //double overlapValue = 81.1;
  double hcoreValue = 0.0e-10;

  if (orbitalA.angularMomentumInt == 0 ) {
    // H_SS ; H_SPx
    if (orbitalB.angularMomentumInt == 0) {
      hcoreValue = parameter_->betaS[orbitalA.element] + \
                   parameter_->betaS[orbitalB.element];
      return hcoreValue * overlapValue * 0.5;
    }else{
      hcoreValue = parameter_->betaS[orbitalA.element] + \
                   parameter_->betaP[orbitalB.element];
      return hcoreValue * overlapValue * 0.5;
    }
  }else{
    //H_PxS ; H_PxPx ; H_PxPy
    if (orbitalB.angularMomentumInt > 0) {
      hcoreValue = parameter_->betaP[orbitalA.element] + \
                   parameter_->betaP[orbitalB.element];
      return hcoreValue * overlapValue * 0.5;
    }else{
      hcoreValue = parameter_->betaP[orbitalA.element] + \
                   parameter_->betaS[orbitalB.element];
      return hcoreValue * overlapValue * 0.5;
    }
  }
}
/***************************************************************************************/ 
double Hcore::ComputeNonDiagonalSameAtom(const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){

  double sumVuvB = 0.0e-10;
  for (size_t i=0;i<infoAOs_->size();++i) {
    if (orbitalA.indexAtom != infoAOs_->orbital[i].indexAtom \
        && infoAOs_->orbital[i].angularMomentumInt == 0) {
      sumVuvB += CoreElectronAttraction(orbitalA,orbitalB,infoAOs_->orbital[i]);
    }
  }
  return sumVuvB;
}
/***************************************************************************************/ 
double Hcore::ComputeDiagonal(const AtomicOrbital& orbitalA){
  if (orbitalA.angularMomentumInt == 0) { // Type orbital s
    return parameter_->uss[orbitalA.element] +\
      ComputeNonDiagonalSameAtom(orbitalA,orbitalA);
  }else{ // Type orbital p
    return parameter_->upp[orbitalA.element] +\
      ComputeNonDiagonalSameAtom(orbitalA,orbitalA);
  }
}
/***************************************************************************************/ 
double Hcore::CoreElectronAttraction(const AtomicOrbital& orbitalAu,\
       const AtomicOrbital& orbitalAv,const AtomicOrbital& orbitalB){
  return -(double)orbitalB.GetCoreCharge() *\
         TwoCenterIntegral::GetValueFromArray(orbitalAu,orbitalAv,orbitalB,orbitalB,\
                                              all2CenterIntegral_);
}
/***************************************************************************************/ 
// OPENACC
/***************************************************************************************/
#ifdef OPENACC_AVIL
void Hcore::ComputeMatrix(){
  cout << "Stop here" << endl;
  cout << "this = " << this  << endl;
  cout << "this.matrixHold_ = " << this->matrixHold_  << endl;
  cout << "this.parameter_ = " << this->parameter_  << endl;
  cout << "this.infoAOs_ = " << this->infoAOs_  << endl;
  cout << "this.all2CenterIntegral_ = " << this->all2CenterIntegral_  << endl;
  cout << "this.overlap_ = " << this->overlap_  << endl;
  cout << "this.overlap_.basisSTO = " << this->overlap_->basisSTO  << endl;

  #pragma acc parallel loop present(this[0:1],infoAOs_,parameter_,all2CenterIntegral_,overlap_,overlap_->basisSTO)
  for (size_t i=0;i<array1DSize_;++i) {
    unsigned int index_ij[2] = {0,0};
    MyMemory::GetIndex_ij_SymetricMatrix(i,index_ij);
    matrixHold_[i] = ComputeElementMatrixLocal(index_ij[0],index_ij[1]);
  }

  Update_hostMatrix();
}
/***************************************************************************************/ 
double Hcore::ComputeElementMatrixLocal(const size_t &i,const size_t &j){
  // Diagonal elements
  if ( i == j ) {
    return ComputeDiagonal(infoAOs_->orbital[i]);
  }
  // Non diagonal element, same atom
  if (infoAOs_->orbital[i].indexAtom == infoAOs_->orbital[j].indexAtom) {
    return ComputeNonDiagonalSameAtom(infoAOs_->orbital[i],infoAOs_->orbital[j]);
  }else{
    return ComputeNonDiagonalDiffAtom(infoAOs_->orbital[i],infoAOs_->orbital[j]);
  }
}
#endif 
/***************************************************************************************/ 
void Hcore::To_device(){
  #pragma acc enter data copyin(this[0:1])
}
/***************************************************************************************/ 
void Hcore::From_device(){
}
/***************************************************************************************/ 
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _HCORE_CPP_
