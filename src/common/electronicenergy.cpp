#ifndef _ELECTRONICENERGY_CPP_
#define _ELECTRONICENERGY_CPP_

#include "mymemory.h"

#include "electronicenergy.h"
/***************************************************************************************/ 
/***************************************************************************************/
ElectronicEnergy::ElectronicEnergy(){
}
/***************************************************************************************/ 
double ElectronicEnergy::ComputeEnergy(const size_t &nAOs,const Hcore &hcore,\
       const FockMatrix &Fmatrix,const DensityMatrix &Pmatrix){
  double energy = 0.0e-10;
  double tmpEnergy;
  // All matrix are symmetric, thus only compute the half lower triangle multiply by two
  // and subtract the diagonal.
  for (size_t i=0;i<nAOs;++i) {
    for (size_t j=0;j<=i;++j) {
      tmpEnergy = 0.0e-10;
      tmpEnergy = hcore.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)] +\
                  Fmatrix.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)];
      tmpEnergy *= Pmatrix.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)];
      energy += tmpEnergy;
    }
  }
  return energy;
}
 
#endif // _ELECTRONICENERGY_CPP_
