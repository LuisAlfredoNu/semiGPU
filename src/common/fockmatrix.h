/*
   Class for compute the fock matrix, the two electron part is computing in situ
*/
#ifndef _FOCKMATRIX_H_
#define _FOCKMATRIX_H_

#include "hcore.h"

#include "basematrix.h"

class FockMatrix : public BaseMatrix {
 public:
  FockMatrix(const Hcore &hcore,const size_t &nAOs);
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods 
  double ComputeElementMatrix(const size_t &i,const size_t &j);
 private:
/***************************************************************************************/ 
  // Variables
  const Hcore* hcore_;
};

 
#endif // _FOCKMATRIX_H_ 
