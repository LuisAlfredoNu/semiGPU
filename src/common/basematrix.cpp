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
  array1DSize_ = MyMemory::GetTotalSizeSymmetricMatrix(Nrows);
}
/***************************************************************************************/ 
BaseMatrix::BaseMatrix(const BaseMatrix &copyof){
  matrixSize_ = copyof.matrixSize_;
  matrixHold_ = copyof.matrixHold_;
  array1DSize_ = copyof.array1DSize_;
}
/***************************************************************************************/ 
bool BaseMatrix::Alloc4Matrix(string matrixName){
  bool isAlloc = MyMemory::AllocSymmetricMatrixReal(matrixName,matrixSize_,matrixHold_);
  To_deviceMatrix();
  return isAlloc;
}
/***************************************************************************************/
bool BaseMatrix::Dealloc4Matrix(){
  From_deviceMatrix();
  bool isDealloc = MyMemory::DeallocSymmetricMatrixReal(matrixHold_);
  return isDealloc;
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
//  size_t total1DArray = matrixSize_ * (matrixSize_ + 1) / 2;
//  unsigned int index_ij[2] = {0,0};
////    std::cout << " Start i , j = " << index_ij[0] << " , " << index_ij[1]  << std::endl;
////    std::cout << " total 1D = " << total1DArray << std::endl;
//#pragma acc parallel loop
//  for (size_t i=0;i<total1DArray;++i) {
//    MyMemory::GetIndex_ij_SymetricMatrix(i,index_ij);
// //   std::cout << "index1D = " << i << " ->  i , j = " << index_ij[0] << " , " << index_ij[1]  << std::endl;
//    matrixHold_[i] = ComputeElementMatrix(index_ij[0],index_ij[1]);
//  }
}
/***************************************************************************************/ 
void BaseMatrix::AssignValue2Matrix(const size_t &i,const size_t &j,const double &value){
  matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)] = value;
}
/***************************************************************************************/ 
void BaseMatrix::To_deviceMatrix(){
  #pragma acc enter data copyin(this[0:1],this->matrixHold_[0:array1DSize_])
  ;
}
void BaseMatrix::From_deviceMatrix(){
  #pragma acc exit data delete(this->matrixHold_[0:array1DSize_],this[0:1])
  ;
}
void BaseMatrix::Update_hostMatrix(){
  #pragma acc update self(this->matrixHold_[0:array1DSize_])
  ;
}
void BaseMatrix::Update_deviceMatrix(){
  #pragma acc update device(this->matrixHold_[0:array1DSize_])
  ;
}
/***************************************************************************************/ 
#endif // _BASEMATRIX_CPP_
