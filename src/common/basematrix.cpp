#ifndef _BASEMATRIX_CPP_
#define _BASEMATRIX_CPP_

#include <string>

#include "mymemory.h"

#include "basematrix.h"
/***************************************************************************************/ 
BaseMatrix::BaseMatrix(){
}
/***************************************************************************************/ 
bool BaseMatrix::Alloc4Matrix(string matrixName,size_t Nrows,double* &matrixHold){
  return MyMemory::AllocSymmetricMatrixReal(matrixName,Nrows,matrixHold);
}
/***************************************************************************************/
bool BaseMatrix::Dealloc4Matrix(double* &matrixHold){
  return MyMemory::DeallocSymmetricMatrixReal(matrixHold);
}
/***************************************************************************************/
void BaseMatrix::ComputeMatrix(double* &matrixHold){}
/***************************************************************************************/ 
void BaseMatrix::AssignValue2Matrix(const size_t &i,const size_t &j,const double &value,\
                                    double* &matrixHold){
  matrixHold[MyMemory::GetIndexSymmetricMatrix(i,j)] = value;
}
#endif // _BASEMATRIX_CPP_
