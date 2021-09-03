#ifndef _SCFCALCULATION_CPP_
#define _SCFCALCULATION_CPP_

#include "mymemory.h"
#include "electronicenergy.h"

#include "scfcalculation.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
SCFCalculation::SCFCalculation(const ListAtomicOrbitals &infoAOs,\
      double**** &all2CenterIntegral,const Hcore &hcore,\
      DensityMatrix &Pmatrix,FockMatrix &Fmatrix){
  infoAOs_ = &infoAOs;
  nAOs_ = infoAOs.orbital.size();
  all2CIntegral_ = all2CenterIntegral;
  hcore_ = &hcore;
  Pmatrix_ = &Pmatrix;
  Fmatrix_ = &Fmatrix;
}
/***************************************************************************************/ 
void SCFCalculation::ComputeSCF(){
  InitialGuess();
  double electronicEnergy = ElectronicEnergy::ComputeEnergy(nAOs_,*hcore_,*Fmatrix_,*Pmatrix_); 
}
/***************************************************************************************/ 
/*
   Plan
   Crear el frirts guess
   Calcular los valores y vectores propios
 */
/***************************************************************************************/ 
void SCFCalculation::InitialGuess(){

  for (size_t i=0;i<nAOs_;++i) {
    Pmatrix_->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,i)] = 1.0;
  }
  Fmatrix_->ComputeMatrix();
  for (size_t i=0;i<nAOs_;++i) {
    Fmatrix_->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,i)] -= 0.5;
  }
}
/***************************************************************************************/ 
void SCFCalculation::GetEigenValVecOfFmatrix(){
}
/***************************************************************************************/ 
#endif // _SCFCALCULATION_CPP_
