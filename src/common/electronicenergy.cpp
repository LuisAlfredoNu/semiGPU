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
#pragma acc parallel loop
  for (size_t i=0;i<nAOs;++i) {
    for (size_t j=0;j<=i;++j) {
      tmpEnergy =   hcore.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)] +\
                  Fmatrix.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)];
      tmpEnergy *= Pmatrix.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)];
      // The other half matrix
      tmpEnergy *= 2.0;
      energy += tmpEnergy;
    }
  }
  // subtract giaonal part
  for (size_t i=0;i<nAOs;++i) {
    tmpEnergy =   hcore.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,i)] +\
                  Fmatrix.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,i)];
    tmpEnergy *= Pmatrix.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,i)];
      energy -= tmpEnergy;
  }
  return 0.5 * energy;
}
 
#endif // _ELECTRONICENERGY_CPP_
