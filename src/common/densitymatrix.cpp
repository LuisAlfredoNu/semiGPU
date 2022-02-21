#ifndef _DENSITYMATRIX_CPP_
#define _DENSITYMATRIX_CPP_

#include "mymemory.h"

#include "densitymatrix.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
DensityMatrix::DensityMatrix(const ListAtomicOrbitals& infoAOs,const double* eigenVec,\
                             const size_t nAOs) : BaseMatrix(nAOs){
  infoAOs_ = &infoAOs;
  eigenVec_ = eigenVec;
  nAOs_ = nAOs;

  nElectron_ = NumberOfElectrons();
  guessSwitch_ = false;
  overlap_ = new Overlap();
}
/***************************************************************************************/ 
double  DensityMatrix::ComputeElementMatrix(const size_t &i,const size_t &j){
  // Compute when the eigenVec are in row-major order 
  // https://en.wikipedia.org/wiki/Row-_and_column-major_order#Column-major_order

  if (guessSwitch_) {
    if (i == j) {
      return (double) infoAOs_->orbital[i].GetCoreCharge() / infoAOs_->orbital[i].GetAOsSize();
    }else{
      double value = 0.0;
      //value = overlap_->ComputeOverlap(infoAOs_->orbital[i],
      //                                 infoAOs_->orbital[j]);
      //value *= (double) infoAOs_->orbital[i].GetCoreCharge();
      //value /= (double) infoAOs_->orbital[i].GetAOsSize() * 2.0;
      return value;
    }
  }else{
    double value = 0.0e-10;
    for (size_t k=0;k<nElectron_/2;++k) {
      // eigenVec_ in row major order
      //value += eigenVec_[nAOs_*i + k] * eigenVec_[nAOs_*j + k];
      // eigenVec_ in column major order
      value += eigenVec_[nAOs_*k + i] * eigenVec_[nAOs_*k + j];
    }
    return value * 2.0;
  }
}
/***************************************************************************************/ 
void DensityMatrix::GuessDensityMatrixSwitch(bool guess){
  guessSwitch_ = guess; 
}
/***************************************************************************************/ 
unsigned int DensityMatrix::NumberOfElectrons(){
  unsigned int nE = 0;
  for (size_t i=0;i<nAOs_;++i) {
    nE += infoAOs_->orbital[i].GetCoreCharge();
    i += infoAOs_->orbital[i].GetAOsSize() - 1  ;
  }

  return nE;
}
/***************************************************************************************/ 
/***************************************************************************************/ 
 
#endif // _DENSITYMATRIX_CPP_
