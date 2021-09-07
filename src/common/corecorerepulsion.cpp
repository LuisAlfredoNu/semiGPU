#ifndef _CORECORERESPULSION_CPP_
#define _CORECORERESPULSION_CPP_

#include <iostream>
using std::cout;
using std::endl;

#include "mymath.h"

#include "corecorerepulsion.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
CoreCoreRepulsion::CoreCoreRepulsion(){
}
/***************************************************************************************/ 
double CoreCoreRepulsion::ComputeRepulsion(const MNDOparameter &parameter,\
                   const ListAtomicOrbitals &infoAOs,double**** all2CenterInt){

  size_t nAOs = infoAOs.orbital.size();
  double energy = 0.0e-10; 
  double tmp_energy;
  for (size_t i=0;i<nAOs;++i) {
    for (size_t j=i;j<nAOs;++j) {
      if (infoAOs.orbital[i].indexAtom != infoAOs.orbital[j].indexAtom ) {
        tmp_energy = (double) (infoAOs.orbital[i].GetCoreCharge() * \
                     infoAOs.orbital[j].GetCoreCharge());
        tmp_energy *= TwoCenterIntegral::GetValueFromArray(
                      infoAOs.orbital[i],infoAOs.orbital[i],
                      infoAOs.orbital[j],infoAOs.orbital[j],all2CenterInt);
        tmp_energy *= AdjustableRepulsion(parameter,infoAOs.orbital[i],infoAOs.orbital[j]);
        energy += tmp_energy;
        //cout << "Hit " << i << " - " << j << "  pair E = " << tmp_energy <<endl;
        if (infoAOs.orbital[j].element > 1) {
          j += 3;
        }
      }
    }
    if (infoAOs.orbital[i].element > 1) {
      i += 3;
    }
  }
  return energy;
}
/***************************************************************************************/ 
double CoreCoreRepulsion::AdjustableRepulsion(const MNDOparameter &parameter,\
       const AtomicOrbital &orbitalA, const AtomicOrbital &orbitalB){
  
  double r_AB = distancePointsV3(orbitalA.coordinates,orbitalB.coordinates);

  if (orbitalB.element == 1 && (orbitalA.element == 7 || orbitalA.element == 8)) {
    return 1.0 + r_AB * exp(-parameter.alpha[orbitalA.element] / convertAngstrom2AU(1.0) *r_AB) + \
           exp(-parameter.alpha[orbitalB.element] / convertAngstrom2AU(1.0) *r_AB) ;
  }
  if (orbitalA.element == 1 && (orbitalB.element == 7 || orbitalB.element == 8)) {
    return 1.0 + r_AB * exp(-parameter.alpha[orbitalB.element] / convertAngstrom2AU(1.0) *r_AB) + \
           exp(-parameter.alpha[orbitalA.element] / convertAngstrom2AU(1.0) *r_AB) ;
  }
    return 1.0 + exp(-parameter.alpha[orbitalA.element] / convertAngstrom2AU(1.0) *r_AB) + \
           exp(-parameter.alpha[orbitalB.element] / convertAngstrom2AU(1.0) *r_AB) ;
}
/***************************************************************************************/ 
#endif // _CORECORERESPULSION_CPP_
