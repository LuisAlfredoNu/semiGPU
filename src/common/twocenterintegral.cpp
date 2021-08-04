
#ifndef _TWOCENTERINTEGRAL_CPP_
#define _TWOCENTERINTEGRAL_CPP_

#include "mymath.h"

#include "twocenterintegral.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

TwoCenterIntegral::TwoCenterIntegral(MNDOparameter& MNDOpara){
  parameter = &MNDOpara;
  multipole = new Multipole(MNDOpara);
}

/***************************************************************************************/ 
double TwoCenterIntegral::ComputeTwoCenterIntegral(const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB,const AtomicOrbital& orbitalC,\
    const AtomicOrbital& orbitalD){

  double integralValue = 0.0e-10;

  // The firts pair of orbitals are in the same atom and the same for the second pair 
  if (orbitalA.nAtom == orbitalB.nAtom && orbitalC.nAtom == orbitalD.nAtom) {
    int pairTypeA = GetPairType(orbitalA,orbitalB);
    int pairTypeB = GetPairType(orbitalC,orbitalD);
    int AOsTypeInt[4] = {orbitalA.angularMomentumInt,orbitalB.angularMomentumInt,\
                         orbitalC.angularMomentumInt,orbitalD.angularMomentumInt};
   
    // If the 4 orbital are in the same atom
    if (orbitalA.nAtom == orbitalC.nAtom) {
      if (pairTypeA == 0) {
        // SS|SS ; SS|PP
        if (pairTypeA == pairTypeB) {
          // SS|SS
          return SelfIntegralTypeSS_SS(orbitalA.element);
        }else{
          // SS|PP
          return SelfIntegralTypeSS_PP(orbitalA.element,AOsTypeInt);
        }
      }else{
        // SP|SP ; PP|PP
        if (pairTypeA == 1) {
          // SP|SP
          return SelfIntegralTypeSP_SP(orbitalA.element,AOsTypeInt);
        }else{
          // PP|PP
          return SelfIntegralTypePP_PP(orbitalA.element,AOsTypeInt);
        }
      }

    } else { // If the two pairs are in differents atoms 
      double atomsDistance = distancePointsV3(orbitalA.coordinates,orbitalC.coordinates);
      double rotationMatrix[3][3];
      if (pairTypeA == 0) {
        // SS|SS ; SS|SP ; SS|PP
        if (pairTypeA == pairTypeB) {
          // SS|SS
          integralValue = IntegralTypeSS_SS(atomsDistance,orbitalA,orbitalC);
        }else{
          // SS|SP ; SS|PP
          if (pairTypeA == 1) {
            // SS|SP
            integralValue = IntegralTypeSS_SP(atomsDistance,orbitalA,orbitalC,\
                                              AOsTypeInt,rotationMatrix);
          }else{
            // SS|PP
          }
        }
      }else{
        // SP|SP ; SP|PP ; PP|PP
        if (pairTypeA == pairTypeB) {
          // SP|SP ; PP|PP
          if (pairTypeA == 1) {
            // SP|SP
          }else{
            // PP|PP
          }
        }else{
          // SP|PP
        }
      }
    }
  }
  return convertHartree2eV(integralValue);
}
/***************************************************************************************/
void TwoCenterIntegral::ComputeRotationMatrix(const double (&vecA)[3],\
    const double (&vecB)[3],double (&matrix)[3][3]){
  double R_ij = distancePointsV3(vecA,vecB);

  for (int i=0;i<3;i++) {
    matrix[2][i] = (vecA[i] - vecB[i]) / R_ij; 
  }

  if (sameReal(matrix[2][0],0.0e-10) && sameReal(matrix[2][1],0.0e-10)) {
    matrix[1][0] = 0.0e-10;
    matrix[1][1] = 1.0;
    matrix[1][2] = 0.0e-10;
    matrix[0][0] = 1.0;
    matrix[0][1] = 0.0e-10;
    matrix[0][2] = 0.0e-10;
  }else{
    double magZ_xy = 1.0 / sqrt(std::pow(matrix[2][0],2.0) + std::pow(matrix[2][1],2.0) );
    matrix[1][0] = magZ_xy * -matrix[2][1];
    matrix[1][1] = magZ_xy * matrix[2][0];
    matrix[1][2] = 0.0e-10;
    crossProductV3(matrix[2],matrix[1],matrix[0]);
  }
}
/***************************************************************************************/
void TwoCenterIntegral::ApplyRotationAOs(int locA, const int& gloA,double& value,\
    const double (&rotMat)[3][3]){

  value = value * rotMat[locA][gloA-1];
}
/***************************************************************************************/ 
double TwoCenterIntegral::IntegralTypeSS_SS(const double& atomsDistance,const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalC){
  return multipole->Interaction_qq(atomsDistance,orbitalA,orbitalC);
}
/***************************************************************************************/ 
double TwoCenterIntegral::IntegralTypeSS_SP(const double& atomsDistance,\
    const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalC,\
    const int (&AOsTypeInt)[4],double (&rotMat)[3][3]){
  
  double value = multipole->Interaction_qq(atomsDistance,orbitalA,orbitalC);
  ComputeRotationMatrix(orbitalA.coordinates,orbitalC.coordinates,rotMat);
  ApplyRotationAOs(0,AOsTypeInt[3],value,rotMat);
  return value;
}
/***************************************************************************************/ 
double TwoCenterIntegral::SelfIntegralTypeSS_SS(const int& element){
  return parameter->gss[element];
}
/***************************************************************************************/ 
double TwoCenterIntegral::SelfIntegralTypeSS_PP(const int& element,const int (&AOsTypeInt)[4]){
  if (AOsTypeInt[2] == AOsTypeInt[3]) {
    return parameter->gsp[element];
  }else{
    return 0.0e-10;
  }
}
/***************************************************************************************/ 
double TwoCenterIntegral::SelfIntegralTypeSP_SP(const int& element,const int (&AOsTypeInt)[4]){
  if (AOsTypeInt[1] == AOsTypeInt[3]) {
    return parameter->hsp[element];
  }else{
    return 0.0e-10;
  }
}
/***************************************************************************************/ 
double TwoCenterIntegral::SelfIntegralTypePP_PP(const int& element,const int (&AOsTypeInt)[4]){
  if (AOsTypeInt[0] == AOsTypeInt[1] && AOsTypeInt[2] == AOsTypeInt[3] ) {
    if (AOsTypeInt[0] == AOsTypeInt[2]) {
      return parameter->gpp[element];
    }else{
      return parameter->gp2[element];
    }
  }
  if (AOsTypeInt[0] == AOsTypeInt[2] && AOsTypeInt[1] == AOsTypeInt[3] ) {
    return 0.5 * (parameter->gpp[element]-parameter->gp2[element]);
  }
  return 0.0e-10;
}
/***************************************************************************************/ 
int TwoCenterIntegral::GetPairType(const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalB){
  return orbitalA.angularMomentum[0] + orbitalA.angularMomentum[1] + orbitalA.angularMomentum[2]\
    + orbitalB.angularMomentum[0] + orbitalB.angularMomentum[1] + orbitalB.angularMomentum[2];
}
/***************************************************************************************/ 
#endif // _TWOCENTERINTEGRAL_H_
