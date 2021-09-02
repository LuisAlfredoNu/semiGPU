#ifndef _FOCKMATRIX_CPP_
#define _FOCKMATRIX_CPP_

#include "mymemory.h"

#include "fockmatrix.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
FockMatrix::FockMatrix(const Hcore &hcore,const size_t &nAOs) : BaseMatrix(nAOs){
  hcore_ = &hcore;
}
/***************************************************************************************/ 
double FockMatrix::ComputeElementMatrix(const size_t &i, const size_t &j){
  return hcore_->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)];
}

#endif // _FOCKMATRIX_CPP_
