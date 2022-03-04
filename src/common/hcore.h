/*
 * Class for compute the Hcore matrix
 * Ref: https://doi.org/10.1021/ja00457a004
 */
#ifndef _HCORE_H_
#define _HCORE_H_

#include "MNDO_parameters.h"
#include "STO-6G.h"
#include "overlap.h"

#include "basematrix.h"

class Hcore : public BaseMatrix {
 public:
  Hcore(const MNDOparameter&,const ListAtomicOrbitals&,double****);
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  double ComputeElementMatrix(const size_t &i,const size_t &j);
#ifdef OPENACC_AVIL
  //OpenACC
  void ComputeMatrix();
  #pragma acc routine seq
  double ComputeElementMatrixLocal(const size_t &i,const size_t &j);
#endif

 private:
/***************************************************************************************/ 
  // Varaibles
  const MNDOparameter* parameter_;
  const ListAtomicOrbitals* infoAOs_;
  double**** all2CenterIntegral_;
  Overlap* overlap_;

/***************************************************************************************/ 
  // Methods
  #pragma acc routine seq
  double ComputeNonDiagonalDiffAtom(const AtomicOrbital&,const AtomicOrbital&);
  #pragma acc routine seq
  double ComputeNonDiagonalSameAtom(const AtomicOrbital&,const AtomicOrbital&);
  #pragma acc routine seq
  double ComputeDiagonal(const AtomicOrbital&);
  #pragma acc routine seq
  double CoreElectronAttraction(const AtomicOrbital&,const AtomicOrbital&,\
         const AtomicOrbital&);
  void To_device();
  void From_device();
};


#endif // _HCORE_H_
