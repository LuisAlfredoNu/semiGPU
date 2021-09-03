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
  SCFCalculation(const ListAtomicOrbitals &,double**** &,const Hcore &,DensityMatrix &,\
                 FockMatrix &);
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  void ComputeSCF();
 private:
/***************************************************************************************/ 
  // Variables
  const ListAtomicOrbitals* infoAOs_;
  size_t nAOs_;
  double**** all2CIntegral_;
  const Hcore* hcore_;
  DensityMatrix* Pmatrix_;
  FockMatrix* Fmatrix_;
/***************************************************************************************/ 
  // Methods
  void InitialGuess();
};

 
#endif // _SCFCALCULATION_H_ 
