#ifndef _HCORE_CPP_
#define _HCORE_CPP_

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
/***************************************************************************************/ 
#endif // _HCORE_CPP_
