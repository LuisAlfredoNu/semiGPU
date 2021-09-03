#ifndef _FOCKMATRIX_CPP_
#define _FOCKMATRIX_CPP_

#include "mymemory.h"
#include "twoelectronmatrix.h"

#include "fockmatrix.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
FockMatrix::FockMatrix(const ListAtomicOrbitals &infoAOs,const Hcore &hcore,\
    const DensityMatrix &Pmatrix,double**** &all2CenterIntegral) : BaseMatrix(infoAOs.orbital.size()){
  infoAOs_ = &infoAOs;
  hcore_ = &hcore;
  Pmatrix_ = &Pmatrix;
  all2CIntegral_ = all2CenterIntegral;
}
/***************************************************************************************/ 
double FockMatrix::ComputeElementMatrix(const size_t &i, const size_t &j){
  double value;
  value = hcore_->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)];
  value += TwoElectronMatrix::ComputeGMatrix_UV(infoAOs_->orbital[i],infoAOs_->orbital[j],\
           *infoAOs_,*Pmatrix_,all2CIntegral_);
  return value;
}

#endif // _FOCKMATRIX_CPP_
