#ifndef _DENSITYMATRIX_CPP_
#define _DENSITYMATRIX_CPP_

#include "mymemory.h"

#include "densitymatrix.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
DensityMatrix::DensityMatrix(const double* eigenVec,const size_t nAOs) : BaseMatrix(nAOs){
  eigenVec_ = eigenVec;
  nAOs_ = nAOs;
}
/***************************************************************************************/ 
double  DensityMatrix::ComputeElementMatrix(const size_t &i,const size_t &j){
  // Compute when the eigenVec are in row-major order 
  // https://en.wikipedia.org/wiki/Row-_and_column-major_order#Column-major_order

  double value = 0.0e-10;
  for (size_t k=0;k<nAOs_/2;++k) {
    // eigenVec_ in row major order
    //value += eigenVec_[nAOs_*i + k] * eigenVec_[nAOs_*j + k];
    // eigenVec_ in column major order
    value += eigenVec_[nAOs_*k + i] * eigenVec_[nAOs_*k + j];
  }
  return value * 2.0;
}
/***************************************************************************************/ 
/***************************************************************************************/ 
 
#endif // _DENSITYMATRIX_CPP_
