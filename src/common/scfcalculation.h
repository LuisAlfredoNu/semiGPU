/*
   Class to compute the SCF procedure 
*/
#ifndef _SCFCALCULATION_H_
#define _SCFCALCULATION_H_

#include "atomicOrbitals.h"
#include "hcore.h"
#include "densitymatrix.h"
#include "fockmatrix.h"

class SCFCalculation {
 public:
  SCFCalculation(const ListAtomicOrbitals &,double**** &,const Hcore &);
/***************************************************************************************/ 
  // Variables
  // Pointer for eigenvectors and eigenvalues
  double* eigenVec;
  double* eigenVal;
  // Pointer to important matrixs
  DensityMatrix* Pmatrix;
  FockMatrix* Fmatrix;
  // Final Energy
  double finalEnergy;
/***************************************************************************************/ 
  // Methods
  bool AllocSCFData();
  bool DeallocSCFData();
  void ComputeSCF();
 private:
/***************************************************************************************/ 
  // Variables
  const ListAtomicOrbitals* infoAOs_;
  size_t nAOs_;
  double**** all2CIntegral_;
  const Hcore* hcore_;
  bool printInfo;
/***************************************************************************************/ 
  // Methods
  void InitialGuess();
  int GetEigenValVecOfFmatrix();
};

 
#endif // _SCFCALCULATION_H_ 
