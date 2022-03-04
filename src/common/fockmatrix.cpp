#ifndef _FOCKMATRIX_CPP_
#define _FOCKMATRIX_CPP_

#include <iostream>
using std::cout;
using std::endl;

#include "mymemory.h"
#include "twoelectronmatrix.h"

#include "fockmatrix.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
FockMatrix::FockMatrix(const ListAtomicOrbitals &infoAOs,const Hcore &hcore,\
    const DensityMatrix &Pmatrix,double**** &all2CenterIntegral) : BaseMatrix(infoAOs.size()){
  infoAOs_ = &infoAOs;
  hcore_ = &hcore;
  Pmatrix_ = &Pmatrix;
  all2CIntegral_ = all2CenterIntegral;

  To_device();
}
/***************************************************************************************/ 
double FockMatrix::ComputeElementMatrix(const size_t &i, const size_t &j){
  double value;
  value = hcore_->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)];
  value += TwoElectronMatrix::ComputeGMatrix_UV(infoAOs_->orbital[i],infoAOs_->orbital[j],\
           *infoAOs_,*Pmatrix_,all2CIntegral_);
  return value;
}
/***************************************************************************************/ 
/***************************************************************************************/
// OPENACC
/***************************************************************************************/
#ifdef OPENACC_AVIL
void FockMatrix::ComputeMatrix(){
  cout << "Stop here" << endl;
  cout << "FockMatrix this = " << this  << endl;
  cout << "FockMatrix this.matrixHold_ = " << this->matrixHold_  << endl;
  cout << "FockMatrix this.infoAOs_ = " << this->infoAOs_  << endl;
  cout << "FockMatrix this.all2CIntegral_ = " << this->all2CIntegral_  << endl;
  cout << "FockMatrix this.hcore_ = " << this->hcore_  << endl;
  cout << "FockMatrix this.Pmatrix_ = " << this->Pmatrix_  << endl;
  #pragma acc parallel loop present(this[0:1],infoAOs_,all2CIntegral_,hcore_,Pmatrix_)
  for (size_t i=0;i<array1DSize_;++i) {
    unsigned int index_ij[2] = {0,0};
    MyMemory::GetIndex_ij_SymetricMatrix(i,index_ij);
    matrixHold_[i] = ComputeElementMatrixLocal(index_ij[0],index_ij[1]);
  }
  Update_hostMatrix();
}
/***************************************************************************************/ 
double FockMatrix::ComputeElementMatrixLocal(const size_t &i, const size_t &j){
  double value;
  size_t index = MyMemory::GetIndexSymmetricMatrix(i,j);
  value = hcore_->matrixHold_[index];
  value += TwoElectronMatrix::ComputeGMatrix_UV(infoAOs_->orbital[i],infoAOs_->orbital[j],\
           *infoAOs_,*Pmatrix_,all2CIntegral_);
  return value;
}
#endif
/***************************************************************************************/ 
void FockMatrix::To_device(){
  #pragma acc enter data copyin(this[0:1])
}
/***************************************************************************************/
void FockMatrix::From_device(){
  #pragma acc exit data delete(this[0:1])
}
/***************************************************************************************/
/***************************************************************************************/
#endif // _FOCKMATRIX_CPP_
