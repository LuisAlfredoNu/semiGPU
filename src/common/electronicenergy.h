/*
   Class for compute the electronic energy from Density matrix, Hcore and Fock matrix
*/
#ifndef _ELECTRONICENERGY_H_
#define _ELECTRONICENERGY_H_

#include "hcore.h"
#include "fockmatrix.h"
#include "densitymatrix.h"

class ElectronicEnergy {
 public:
  ElectronicEnergy();
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  static double ComputeEnergy(const size_t &nAOs,const Hcore &,const FockMatrix &,\
         const DensityMatrix &);
 private:
};

#endif // _ELECTRONICENERGY_H_ 
