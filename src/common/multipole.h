/*
 * Class for compute the multipole interactions  
 * Ref  https://doi.org/10.1007/BF00548085
 */

#ifndef _MULTIPOLE_H_
#define _MULTIPOLE_H_

#include "atomicOrbitals.h"
#include "MNDO_parameters.h"

/***************************************************************************************/ 
class Multipole{
/***************************************************************************************/ 
 public:
  Multipole(MNDOparameter&);
/***************************************************************************************/ 
  // Variable
  MNDOparameter* parameter;
/***************************************************************************************/ 
  // Methods
  double Interaction_qq(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_qUz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_Uzq(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_qQpipi(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_Qpipiq(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_qQzz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_Qzzq(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_UpiUpi(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_UzUz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_UpiQpiz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QpizUpi(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_UzQpipi(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QpipiUz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_UzQz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QzUz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QpipiQpipi(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QxxQyy(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QpipiQzz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QzzQpipi(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QzzQzz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QpizQpiz(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);
  double Interaction_QxyQxy(const double& atomDistance,const AtomicOrbital&, const AtomicOrbital&);

 private:
  void AditiveTerms(double&,const int&,int,const int&,int);
};
#endif // _MULTIPOLE_H_
