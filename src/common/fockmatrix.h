/*
   Class for compute the fock matrix, the two electron part is computing in situ
*/
#ifndef _FOCKMATRIX_H_
#define _FOCKMATRIX_H_

#include "atomicOrbitals.h"
#include "hcore.h"
#include "densitymatrix.h"
#include "twocenterintegral.h"

#include "basematrix.h"

class FockMatrix : public BaseMatrix {
 public:
  FockMatrix(const ListAtomicOrbitals &,const Hcore &, const DensityMatrix &, double**** &);
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods 
  double ComputeElementMatrix(const size_t &i,const size_t &j);
#ifdef OPENACC_AVIL
  //OpenACC
  #pragma acc routine 
  double ComputeElementMatrixLocal(const size_t &i,const size_t &j);
  void ComputeMatrix();
#endif // OPENACC_AVIL
 private:
/***************************************************************************************/ 
  // Variables
  const ListAtomicOrbitals* infoAOs_;
  const Hcore* hcore_;
  const DensityMatrix* Pmatrix_;
  double**** all2CIntegral_;
/***************************************************************************************/ 
  // Methods
  //OpenACC
  void To_device();
  void From_device();
};

 
#endif // _FOCKMATRIX_H_ 
