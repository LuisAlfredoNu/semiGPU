#ifndef _OVERLAP_CPP_
#define _OVERLAP_CPP_


#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <cmath>

#include "atomicOrbitals.h"
#include "mymemory.h"
#include "mymath.h"

#include "overlap.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

Overlap::Overlap(){
}

bool Overlap::GetOverlapMatrix(const vector<AtomicOrbital> infoOrbitals,double* &overlapMatrix){
  
  string ovMtx = "overlapMatrix";
  int NOrbitals = infoOrbitals.size();

  cout << "overlapMatrix " << overlapMatrix << endl;
  bool correctAlloc = MyMemory::AllocSymetricMatrixReal(ovMtx,NOrbitals,overlapMatrix);
  cout << "overlapMatrix " << overlapMatrix << endl;

  if (! correctAlloc) {
    return false;
  }

  double overlapValue = 0.0;
  
  for (int i=0;i<NOrbitals;i++) {
    for (int j=0;j <= i;j++) {
      // TODO 
      // - Checar si regresar un array 1D o convertir AllocSymetricMatrixReal a un array 2D
      overlapValue = ComputeOverlap(infoOrbitals[i],infoOrbitals[j]);
      overlapMatrix[MyMemory::GetIndexSymetrixMatrix(i,j)] = overlapValue;
    }
  }

  return true;
}
/***************************************************************************************/ 
double Overlap::ComputeOverlap(const AtomicOrbital orbitalA,const AtomicOrbital orbitalB){

  double overlapGTO = 0.0;
  double overlapSTO = 0.0e-10;

  int atomTypeA = orbitalA.element;
  int atomTypeB = orbitalB.element;

  double atomDistance = distancePointsV3(orbitalA.coordinates,orbitalB.coordinates);

  // Multiply and sum of exponents
  double a_times_a = 0.0;
  double a_plus_a = 0.0;

  for(int gtoA = 0; gtoA < 6;gtoA++){
    for(int gtoB = 0; gtoB < 6;gtoB++){
      a_times_a = basisSTO.value[atomTypeA].exponent[gtoA] * \
                  basisSTO.value[atomTypeB].exponent[gtoB];
      a_plus_a = basisSTO.value[atomTypeA].exponent[gtoA] +  \
                 basisSTO.value[atomTypeB].exponent[gtoB];

      overlapGTO = 2.0 * sqrt(a_times_a)/a_plus_a;
      overlapGTO = sqrt(overlapGTO);
      overlapGTO = overlapGTO * overlapGTO * overlapGTO;
      overlapGTO *= exp(-a_times_a / a_plus_a * atomDistance * atomDistance);

      // normalization factor

      overlapSTO += overlapGTO * basisSTO.value[atomTypeA].coeffS[gtoA] * \
                    basisSTO.value[atomTypeB].coeffS[gtoB];
    }
  }
  return overlapSTO;
}


#endif // _OVERLAP_CPP_
