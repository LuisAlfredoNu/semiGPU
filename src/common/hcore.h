/*
 * Class for compute the Hcore matrix
 * Ref: https://doi.org/10.1021/ja00457a004
 */
#ifndef _HCORE_H_
#define _HCORE_H_

#include "MNDO_parameters.h"
#include "overlap.h"

#include "basematrix.h"

class Hcore : public BaseMatrix {
 public:
  Hcore(const MNDOparameter&,const ListAtomicOrbitals&,double****);
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  void ComputeMatrix(double* &hcoreMatrix);

 private:
/***************************************************************************************/ 
  // Varaibles
  const MNDOparameter* parameter;
  const ListAtomicOrbitals* infoAOs;
  double**** all2CenterIntegral;
  Overlap* overlap;

/***************************************************************************************/ 
  // Methods
  double ComputeNonDiagonalDiffAtom(const AtomicOrbital&,const AtomicOrbital&);
  double ComputeNonDiagonalSameAtom(const AtomicOrbital&,const AtomicOrbital&);
  double ComputeDiagonal(const AtomicOrbital&);
  double CoreElectronAttraction(const AtomicOrbital&,const AtomicOrbital&,\
         const AtomicOrbital&);
};


#endif // _HCORE_H_
