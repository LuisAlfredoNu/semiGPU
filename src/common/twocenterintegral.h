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
/*********************W******************************************************************/ 
  // Methods
  double ComputeTwoCenterIntegral(const AtomicOrbital& orbitalA,\
      const AtomicOrbital& orbitalB,const AtomicOrbital& orbitalC,\
      const AtomicOrbital& orbitalD);
  // All data for 2 center integral
  static bool Alloc4AllTwoCenterIntegral(const vector<Atom>& molecule,\
      double**** &all2CenterIntegral);
  static bool Dealloc4AllTwoCenterIntegral(const vector<Atom>& molecule,\
      double**** &all2CenterIntegral);
  void ComputeAllTwoCenterIntegral(const ListAtomicOrbitals& infoAOs,\
      double**** &all2CenterIntegral);
  static double GetValueFromArray(const AtomicOrbital&,const AtomicOrbital&,const AtomicOrbital&,const AtomicOrbital&,double**** all2CenterIntegral);
/***************************************************************************************/
/***************************************************************************************/ 
 private:
/***************************************************************************************/ 
  // Varaibles
  MNDOparameter* parameter;
  Multipole* multipole;
/***************************************************************************************/ 
  // Methods
  // Diff Atom
  double IntegralTypeSS_SS(const double&,const AtomicOrbital&,const AtomicOrbital&);
  double IntegralTypeSS_SP(const double&,const AtomicOrbital&,const AtomicOrbital&,\
      const int (&AOsTypeInt)[4],double (&rotMat)[3][3]);
  double IntegralTypeSS_PP(const double&,const AtomicOrbital&,const AtomicOrbital&,\
      const int (&AOsTypeInt)[4],double (&rotMat)[3][3]);
  double IntegralTypeSP_SP(const double&,const AtomicOrbital&,const AtomicOrbital&,\
      const int (&AOsTypeInt)[4],double (&rotMat)[3][3]);
  double IntegralTypeSP_PP(const double&,const AtomicOrbital&,const AtomicOrbital&,\
      const int (&AOsTypeInt)[4],double (&rotMat)[3][3]);
  double IntegralTypePP_PP(const double&,const AtomicOrbital&,const AtomicOrbital&,\
      const int (&AOsTypeInt)[4],double (&rotMat)[3][3]);
  // Same atoms
  double SelfIntegralTypeSS_SS(const int&);
  double SelfIntegralTypeSS_PP(const int&, const int (&AOsTypeInt)[4]);
  double SelfIntegralTypePP_SS(const int&, const int (&AOsTypeInt)[4]);
  double SelfIntegralTypeSP_SP(const int&, const int (&AOsTypeInt)[4]);
  double SelfIntegralTypePP_PP(const int&, const int (&AOsTypeInt)[4]);
  // Functions for type integral
  int GetPairType(const AtomicOrbital&,const AtomicOrbital&);
  double CorrectOrderIntegral(int&,int&,int (&)[4]);

  // Rotations
  void ComputeRotationMatrix(const double (&vecA)[3],const double (&vecB)[3],\
      double (&matrix)[3][3]);
  double ApplyRotationAOs(int,const int&,const double&, const double (&rotMat)[3][3]);
  double ApplyRotationAOs(int,int,const int&,const int&,const double&,const double (&rotMat)[3][3]);
  double ApplyRotationAOs(int,int,int,const int&,const int&,const int&,const double&,const double (&rotMat)[3][3]);
  double ApplyRotationAOs(int,int,int,int,const int&,const int&,const int&,const int&,const double&,const double (&rotMat)[3][3]);

  // Funtions fro compute by atoms pairs
  void ComputePair_HH(const int&,int&,const vector<AtomicOrbital>&,double**&);
  void ComputePair_HX(      int&,int&,const vector<AtomicOrbital>&,double**&);
  void ComputePair_XX(const int&,int&,const vector<AtomicOrbital>&,double**&);
};
#endif // _TWOCENTERINTEGRAL_H_
