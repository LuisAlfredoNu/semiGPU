/*
   Class for compute the core-core repultion in all system
*/
#ifndef _CORECORERESPULSION_H_
#define _CORECORERESPULSION_H_

#include "atomicOrbitals.h"
#include "twocenterintegral.h"
#include "MNDO_parameters.h"

class CoreCoreRepulsion {
 public:
  CoreCoreRepulsion();
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  static double ComputeRepulsion(const MNDOparameter &,const ListAtomicOrbitals &,
         double**** );
 private:
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  static double AdjustableRepulsion(const MNDOparameter &,const AtomicOrbital &,\
         const AtomicOrbital &);

};

 
#endif // _CORECORERESPULSION_H_ 
