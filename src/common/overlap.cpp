#ifndef _OVERLAP_CPP_
#define _OVERLAP_CPP_


#include <string>
using std::string;
#include <cmath>

#include "atomicOrbitals.h"
#include "mymemory.h"
#include "mymath.h"
#include "STO-6G.h"

#include "overlap.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
Overlap::Overlap() : BaseMatrix(){
  basisSTO = new STO_6G();
}
/***************************************************************************************/ 
Overlap::Overlap(ListAtomicOrbitals infoAOs) : BaseMatrix(infoAOs.orbital.size()){
  basisSTO = new STO_6G();
  infoAOs_ = &infoAOs;
}
/***************************************************************************************/ 
double Overlap::ComputeElementMatrix(const size_t &i,const size_t &j){
  return ComputeOverlap(infoAOs_->orbital[i],infoAOs_->orbital[j]);
}
/***************************************************************************************/ 
double Overlap::ComputeOverlap(const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalB){
  int atomTypeA = orbitalA.element;
  int atomTypeB = orbitalB.element;

  int sumAngularMomentumA = orbitalA.angularMomentum[0] + orbitalA.angularMomentum[1] + \
                            orbitalA.angularMomentum[2];
  int sumAngularMomentumB = orbitalB.angularMomentum[0] + orbitalB.angularMomentum[1] + \
                            orbitalB.angularMomentum[2];

  int xyzAngularMomentumA= 0, xyzAngularMomentumB = 0;

  if (sumAngularMomentumA == 1 ) {
    if (orbitalA.angularMomentum[0] == 1) {
      xyzAngularMomentumA = 0;
    } else if (orbitalA.angularMomentum[1] == 1){
      xyzAngularMomentumA = 1;
    } else if (orbitalA.angularMomentum[2] == 1) {
      xyzAngularMomentumA = 2;
    }
  }
  if (sumAngularMomentumB == 1 ) {
    if (orbitalB.angularMomentum[0] == 1) {
      xyzAngularMomentumB = 0;
    } else if (orbitalB.angularMomentum[1] == 1){
      xyzAngularMomentumB = 1;
    } else if (orbitalB.angularMomentum[2] == 1) {
      xyzAngularMomentumB = 2;
    }
  }

  double atomDistance = distancePointsV3(orbitalA.coordinates,orbitalB.coordinates);
  double rABx = 0.0, rABy = 0.0;

  rABx = orbitalA.coordinates[xyzAngularMomentumA] - orbitalB.coordinates[xyzAngularMomentumA] ;
  rABy = orbitalA.coordinates[xyzAngularMomentumB] - orbitalB.coordinates[xyzAngularMomentumB] ;

  double overlapSTO = 0.0e-10;
  double overlapGTO = 0.0;
  // Multiply and sum of exponents
  double a_times_b = 0.0;
  double a_plus_b = 0.0;
  // Normalization factor
  double normaFactor;

  for(int gtoA = 0; gtoA < 6;gtoA++){
    for(int gtoB = 0; gtoB < 6;gtoB++){
      a_times_b = basisSTO->value[atomTypeA].exponent[gtoA] * \
                  basisSTO->value[atomTypeB].exponent[gtoB];
      a_plus_b = basisSTO->value[atomTypeA].exponent[gtoA] +  \
                 basisSTO->value[atomTypeB].exponent[gtoB];

      Overlap_SS(overlapGTO,a_times_b,a_plus_b,atomDistance);

      if (sumAngularMomentumA == 1 && sumAngularMomentumB == 0) {
        Overlap_PS(overlapGTO,basisSTO->value[atomTypeA].exponent[gtoA],\
                              basisSTO->value[atomTypeB].exponent[gtoB],\
                              a_plus_b,rABx);
      
        normaFactor = basisSTO->value[atomTypeA].coeffP[gtoA] * \
                      basisSTO->value[atomTypeB].coeffS[gtoB];
      }else if (sumAngularMomentumA == 0 && sumAngularMomentumB == 1) {
        Overlap_SP(overlapGTO,basisSTO->value[atomTypeA].exponent[gtoA],\
                              basisSTO->value[atomTypeB].exponent[gtoB],\
                              a_plus_b,rABy);

        normaFactor = basisSTO->value[atomTypeA].coeffS[gtoA] * \
                      basisSTO->value[atomTypeB].coeffP[gtoB];
      }else if (sumAngularMomentumA == 1 && sumAngularMomentumB == 1) {
        if (xyzAngularMomentumA == xyzAngularMomentumB) {
          Overlap_PxPx(overlapGTO,a_times_b,a_plus_b,rABx);
        } else{
          Overlap_PxPy(overlapGTO,a_times_b,a_plus_b,rABx,rABy);
        }
        normaFactor = basisSTO->value[atomTypeA].coeffP[gtoA] * \
                      basisSTO->value[atomTypeB].coeffP[gtoB];
      } else{
        normaFactor = basisSTO->value[atomTypeA].coeffS[gtoA] * \
                      basisSTO->value[atomTypeB].coeffS[gtoB];
      }
      overlapSTO += normaFactor * overlapGTO;
    }
  }
  return overlapSTO;
}
/***************************************************************************************/ 
void Overlap::Overlap_SS(double& overlapGTO,const double& a_times_b,const double& a_plus_b,\
     const double& atomDistance){
  overlapGTO = 2.0 * sqrt(a_times_b)/a_plus_b;
  overlapGTO = sqrt(overlapGTO);
  overlapGTO = overlapGTO * overlapGTO * overlapGTO;
  overlapGTO *= exp(-a_times_b / a_plus_b * atomDistance * atomDistance);
}
/***************************************************************************************/ 
void Overlap::Overlap_PS(double& overlapGTO,const double& expA,const double& expB,const double& a_plus_b,const double& rAB){
  overlapGTO *= - 2.0 * sqrt(expA) * expB / a_plus_b;
  overlapGTO *= rAB;
}
/***************************************************************************************/ 
void Overlap::Overlap_SP(double& overlapGTO,const double& expA,const double& expB,const double& a_plus_b,const double& rAB){
  overlapGTO *= 2.0 * expA * sqrt(expB) / a_plus_b;
  overlapGTO *= rAB;
}
/***************************************************************************************/ 
void Overlap::Overlap_PxPy(double& overlapGTO,const double& a_times_b, const double& a_plus_b, \
     const double& rABx, const double& rABy){
  overlapGTO *= - 4.0 * sqrt(a_times_b * a_times_b * a_times_b) / (a_plus_b * a_plus_b);
  overlapGTO *= rABx * rABy ;
}
/***************************************************************************************/ 
void Overlap::Overlap_PxPx(double& overlapGTO,const double& a_times_b, const double& a_plus_b, \
     const double& rAB){
  overlapGTO *= 4.0 * sqrt(a_times_b) / (a_plus_b);
  overlapGTO *= 0.5 - rAB * rAB * a_times_b / a_plus_b;
}

#endif // _OVERLAP_CPP_
