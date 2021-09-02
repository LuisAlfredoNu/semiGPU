#ifndef _DENSITYMATRIX_H_
#define _DENSITYMATRIX_H_

#include "basematrix.h"

class DensityMatrix : public BaseMatrix {
 public:
  DensityMatrix(const double* eigenVec,const size_t nAOs);
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  double ComputeElementMatrix(const size_t &i,const size_t &j);
 private:
/***************************************************************************************/ 
  // Variables
  size_t nAOs_;
  const double* eigenVec_;

/***************************************************************************************/ 
  // Methods
};

 
#endif // _DENSITYMATRIX_H_ 
