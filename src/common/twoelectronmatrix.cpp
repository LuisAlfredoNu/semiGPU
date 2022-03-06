#ifndef _TWOALECTRONMATRIX_CPP_
#define _TWOALECTRONMATRIX_CPP_

#include "mymath.h"
#include "mymemory.h"

#include "twoelectronmatrix.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
TwoElectronMatrix::TwoElectronMatrix(){
}
/***************************************************************************************/ 
/***************************************************************************************/ 
double TwoElectronMatrix::ComputeGMatrix_UV(const AtomicOrbital &orbitalU,\
       const AtomicOrbital &orbitalV,const ListAtomicOrbitals &infoAOs,\
       const DensityMatrix &Pmatrix,double**** &all2CenterIntegral){

  double value = 0.0e-10;
  double valueTmp = 0.0e-10;
  #pragma acc loop independent private(valueTmp,value)  reduction(+:value) 
  for (size_t i=0;i<infoAOs.size();++i) {
    #pragma acc loop independent private(valueTmp) vector
    for (size_t j=0;j<infoAOs.size();++j) {
      valueTmp = 0.0e-10;
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
