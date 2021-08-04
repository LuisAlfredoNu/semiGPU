/* 
 * Class for compute the integral for two center on MNDO model
 * Ref: https://doi.org/10.1007/BF00548085 
 */

#ifndef _TWOCENTERINTEGRAL_H_
#define _TWOCENTERINTEGRAL_H_

#include "atomicOrbitals.h"
#include "MNDO_parameters.h"
#include "multipole.h"

class TwoCenterIntegral{
 public:
  TwoCenterIntegral(MNDOparameter&);
/***************************************************************************************/ 
  // Varaibles
  MNDOparameter* parameter;
  Multipole* multipole;
/*********************W******************************************************************/ 
  // Methods
  double ComputeTwoCenterIntegral(const AtomicOrbital& orbitalA,\
      const AtomicOrbital& orbitalB,const AtomicOrbital& orbitalC,\
      const AtomicOrbital& orbitalD);
 private:
/***************************************************************************************/ 
  // Varaibles
/***************************************************************************************/ 
  // Methods
  // Diff Atom
  double IntegralTypeSS_SS(const double&,const AtomicOrbital&,const AtomicOrbital&);
  double IntegralTypeSS_SP(const double&,const AtomicOrbital&,const AtomicOrbital&,\
      const int (&AOsTypeInt)[4],double (&rotMat)[3][3]);
  // Same atoms
  double SelfIntegralTypeSS_SS(const int&);
  double SelfIntegralTypeSS_PP(const int&, const int (&AOsTypeInt)[4]);
  double SelfIntegralTypeSP_SP(const int&, const int (&AOsTypeInt)[4]);
  double SelfIntegralTypePP_PP(const int&, const int (&AOsTypeInt)[4]);
  int GetPairType(const AtomicOrbital&,const AtomicOrbital&);

  // Rotations
  void ComputeRotationMatrix(const double (&vecA)[3],const double (&vecB)[3],\
      double (&matrix)[3][3]);
  void ApplyRotationAOs(int,const int&,double&, const double (&rotMat)[3][3]);

};
#endif // _TWOCENTERINTEGRAL_H_