#ifndef _HCORE_CPP_
#define _HCORE_CPP_

#include "atomicOrbitals.h"
#include "twocenterintegral.h"
#include "mymemory.h"

#include "hcore.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
Hcore::Hcore(const MNDOparameter& MNDOpara,const ListAtomicOrbitals& allAOs,\
             double**** all2CenInt){
  parameter = &MNDOpara;
  infoAOs = &allAOs;
  all2CenterIntegral = all2CenInt;
  overlap = new Overlap();
}
/***************************************************************************************/ 
void Hcore::ComputeMatrix(double* &hcoreMatrix){
  for (size_t i=0;i<infoAOs->orbital.size();++i) {
    for (size_t j=0;j <= i ;++j) {
      // Diagonal elements
      if ( i == j ) {
          AssignValue2Matrix(i,j,ComputeDiagonal(infoAOs->orbital[i]),hcoreMatrix);
        continue;
      }
      // Non diagonal element, same atom
      if (infoAOs->orbital[i].indexAtom == infoAOs->orbital[j].indexAtom) {
        AssignValue2Matrix(i,j,\
                  ComputeNonDiagonalSameAtom(infoAOs->orbital[i],infoAOs->orbital[j]),\
                  hcoreMatrix);
      }else{
        AssignValue2Matrix(i,j,\
                  ComputeNonDiagonalDiffAtom(infoAOs->orbital[i],infoAOs->orbital[j]),\
                  hcoreMatrix);
      }
    }
  }
}
/***************************************************************************************/ 
double Hcore::ComputeNonDiagonalDiffAtom(const AtomicOrbital& orbitalA,\
       const AtomicOrbital& orbitalB){

  double overlapValue = overlap->ComputeOverlap(orbitalA,orbitalB);
  double hcoreValue = 0.0e-10;


  if (orbitalA.angularMomentumInt == 0 ) {
    // H_SS ; H_SPx
    if (orbitalB.angularMomentumInt == 0) {
      hcoreValue = parameter->betaS[orbitalA.element] + \
                   parameter->betaS[orbitalB.element];
      return hcoreValue * overlapValue * 0.5;
    }else{
      hcoreValue = parameter->betaS[orbitalA.element] + \
                   parameter->betaP[orbitalB.element];
      return hcoreValue * overlapValue * 0.5;
    }
  }else{
    //H_PxS ; H_PxPx ; H_PxPy
    if (orbitalB.angularMomentumInt > 0) {
      hcoreValue = parameter->betaP[orbitalA.element] + \
                   parameter->betaP[orbitalB.element];
      return hcoreValue * overlapValue * 0.5;
    }else{
      hcoreValue = parameter->betaP[orbitalA.element] + \
                   parameter->betaS[orbitalB.element];
      return hcoreValue * overlapValue * 0.5;
    }
  }
}
/***************************************************************************************/ 
double Hcore::ComputeNonDiagonalSameAtom(const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){

  double sumVuvB = 0.0e-10;
  for (size_t i=0;i<infoAOs->orbital.size();++i) {
    if (orbitalA.indexAtom != infoAOs->orbital[i].indexAtom \
        && infoAOs->orbital[i].angularMomentumInt == 0) {
      sumVuvB += CoreElectronAttraction(orbitalA,orbitalB,infoAOs->orbital[i]);
    }
  }
  return sumVuvB;
}
/***************************************************************************************/ 
double Hcore::ComputeDiagonal(const AtomicOrbital& orbitalA){
  if (orbitalA.angularMomentumInt == 0) { // Type orbital s
    return parameter->uss[orbitalA.element] +\
      ComputeNonDiagonalSameAtom(orbitalA,orbitalA);
  }else{ // Type orbital p
    return parameter->upp[orbitalA.element] +\
      ComputeNonDiagonalSameAtom(orbitalA,orbitalA);
  }
}
/***************************************************************************************/ 
double Hcore::CoreElectronAttraction(const AtomicOrbital& orbitalAu,\
       const AtomicOrbital& orbitalAv,const AtomicOrbital& orbitalB){
  return -(double)orbitalB.GetCoreCharge() *\
         TwoCenterIntegral::GetValueFromArray(orbitalAu,orbitalAv,orbitalB,orbitalB,\
                                              all2CenterIntegral);
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _HCORE_CPP_
