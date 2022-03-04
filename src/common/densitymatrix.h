#ifndef _DENSITYMATRIX_H_
#define _DENSITYMATRIX_H_

#include "overlap.h"
#include "atomicOrbitals.h"

#include "basematrix.h"

class DensityMatrix : public BaseMatrix {
 public:
  DensityMatrix(const ListAtomicOrbitals& infoAOs,const double* eigenVec,const size_t nAOs);
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  double ComputeElementMatrix(const size_t &i,const size_t &j);
  void GuessDensityMatrixSwitch(bool);
#ifdef OPENACC_AVIL
  //OpenACC
  void ComputeMatrix();
  double ComputeElementMatrixLocal(const size_t &i,const size_t &j);
#endif
 private:
/***************************************************************************************/ 
  // Variables
  size_t nAOs_;
  unsigned int nElectron_;
  const double* eigenVec_;
  bool guessSwitch_;
  const ListAtomicOrbitals* infoAOs_;
  //Overlap* overlap_;

/***************************************************************************************/ 
  // Methods
  unsigned int NumberOfElectrons();
  void To_device();
  void From_device();
};

 
#endif // _DENSITYMATRIX_H_ 
