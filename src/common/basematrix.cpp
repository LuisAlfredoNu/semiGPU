#ifndef _BASEMATRIX_CPP_
#define _BASEMATRIX_CPP_

#include "mymemory.h"

#include "basematrix.h"
/***************************************************************************************/ 
BaseMatrix::BaseMatrix(){
}
/***************************************************************************************/ 
BaseMatrix::BaseMatrix(size_t Nrows){
  matrixSize_ = Nrows;
  matrixHold_ = NULL;
}
/***************************************************************************************/ 
bool BaseMatrix::Alloc4Matrix(string matrixName){
  return MyMemory::AllocSymmetricMatrixReal(matrixName,matrixSize_,matrixHold_);
}
/***************************************************************************************/
bool BaseMatrix::Dealloc4Matrix(){
  return MyMemory::DeallocSymmetricMatrixReal(matrixHold_);
}
/***************************************************************************************/
double BaseMatrix::ComputeElementMatrix(const size_t &i,const size_t &j){
  return -1.0;
}
/***************************************************************************************/ 
void BaseMatrix::ComputeMatrix(){
  for (size_t i=0;i<matrixSize_;++i) {
    for (size_t j=0;j<=i;++j) {
      AssignValue2Matrix(i,j,ComputeElementMatrix(i,j));
    }
  }
}
/***************************************************************************************/ 
void BaseMatrix::AssignValue2Matrix(const size_t &i,const size_t &j,const double &value){
  matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)] = value;
}
#endif // _BASEMATRIX_CPP_
