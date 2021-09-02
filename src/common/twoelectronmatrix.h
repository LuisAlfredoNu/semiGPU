/*
   Compute two electron matrix or G matrix
   Ref:  arXiv:1806.06147v3  
         Ec: 25
*/
#ifndef _TWOALECTRONMATRIX_H_
#define _TWOALECTRONMATRIX_H_

#include "atomicOrbitals.h"
#include "densitymatrix.h"

class TwoElectronMatrix {
 public:
  TwoElectronMatrix();
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  static double ComputeGMatrix_UV(const AtomicOrbital &orbitalU,const AtomicOrbital &orbitalV,\
         const ListAtomicOrbitals &infoAOs,double**** &all2CenterIntegral,\
         const DensityMatrix &Pmatrix);
 private:
};

 
#endif // _TWOALECTRONMATRIX_H_ 