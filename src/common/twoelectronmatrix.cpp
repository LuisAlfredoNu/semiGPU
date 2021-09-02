#ifndef _TWOALECTRONMATRIX_CPP_
#define _TWOALECTRONMATRIX_CPP_

#include "mymath.h"
#include "mymemory.h"
#include "twocenterintegral.h"

#include "twoelectronmatrix.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
TwoElectronMatrix::TwoElectronMatrix(){
}
/***************************************************************************************/ 
/***************************************************************************************/ 
double TwoElectronMatrix::ComputeGMatrix_UV(const AtomicOrbital &orbitalU,\
       const AtomicOrbital &orbitalV,const ListAtomicOrbitals &infoAOs,\
       double**** &all2CenterIntegral,const DensityMatrix &Pmatrix){
  double value = 0.0e-10;
  double valueTmp = 0.0e-10;

  for (size_t i=0;i<infoAOs.orbital.size();++i) {
    for (size_t j=0;j<infoAOs.orbital.size();++j) {
      valueTmp = TwoCenterIntegral::GetValueFromArray(orbitalU,orbitalV,\
              infoAOs.orbital[i],infoAOs.orbital[j],all2CenterIntegral);
      valueTmp -= 0.5 * TwoCenterIntegral::GetValueFromArray(orbitalU,infoAOs.orbital[i],\
               orbitalV,infoAOs.orbital[j],all2CenterIntegral);
      valueTmp *= Pmatrix.matrixHold_[MyMemory::GetIndexFullSymmetricMatrix(i,j)];

      value += valueTmp;
    }
  }
  return value;
}
/***************************************************************************************/ 
/***************************************************************************************/ 

#endif // _TWOALECTRONMATRIX_CPP_
