
#ifndef _TWOCENTERINTEGRAL_CPP_
#define _TWOCENTERINTEGRAL_CPP_

#include "iostream"
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "mymath.h"
#include "mymemory.h"
#include "atom.h"

#include "twocenterintegral.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

TwoCenterIntegral::TwoCenterIntegral(MNDOparameter& MNDOpara){
  parameter = &MNDOpara;
  multipole = new Multipole(MNDOpara);
}
/***************************************************************************************/ 
bool TwoCenterIntegral::Alloc4AllTwoCenterIntegral(const vector<Atom>& molecule,\
    double**** &all2CenterIntegral){

  int NAtoms = molecule.size();
  // Alloc for each atom
  if (!(all2CenterIntegral = new double*** [NAtoms])){
    return false;
  }
  // Alloc for the second atom
  for (int i=1;i <= NAtoms;i++) {
    if (!(all2CenterIntegral[i-1] = new double** [i])){
      return false;
    }
  }
  bool dataStatus = false;
  for (int i=0;i<NAtoms;i++) {
    for (int j=0;j <= i;j++) {
      if (molecule[i].atomNumber == 1 && molecule[j].atomNumber  == 1) {
        // Alloc for H-H pair
        dataStatus = MyMemory::Alloc2DRealArray("all2CenterIntegral",1,1,all2CenterIntegral[i][j],0.0);
        continue;
      }
      if (molecule[i].atomNumber > 1 && molecule[j].atomNumber > 1) {
        // Alloc for X-X pair
        dataStatus = MyMemory::Alloc2DRealArray("all2CenterIntegral",10,10,all2CenterIntegral[i][j],0.0);
        continue;
      }else{
        // Alloc for H-X pair
        dataStatus = MyMemory::Alloc2DRealArray("all2CenterIntegral",1,10,all2CenterIntegral[i][j],0.0);
        continue;
      }
    }
  }
  return dataStatus;
}
/***************************************************************************************/ 
double TwoCenterIntegral::ComputeTwoCenterIntegral(const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB,const AtomicOrbital& orbitalC,\
    const AtomicOrbital& orbitalD){

  double integralValue = 0.0e-10;

  // The firts pair of orbitals are in the same atom and the same for the second pair 
  if (orbitalA.indexAtom == orbitalB.indexAtom && orbitalC.indexAtom == orbitalD.indexAtom) {
    int pairTypeA = GetPairType(orbitalA,orbitalB);
    int pairTypeB = GetPairType(orbitalC,orbitalD);
    int AOsTypeInt[4] = {orbitalA.angularMomentumInt,orbitalB.angularMomentumInt,\
                         orbitalC.angularMomentumInt,orbitalD.angularMomentumInt};
  
    // If the 4 orbital are in the same atom
    if (orbitalA.indexAtom == orbitalC.indexAtom) {
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
      double invertIntegral = CorrectOrderIntegral(pairTypeA,pairTypeB,AOsTypeInt);
      if (pairTypeA == 0) {
        // SS|SS ; SS|SP ; SS|PP
        if (pairTypeA == pairTypeB) {
          // SS|SS
          integralValue = IntegralTypeSS_SS(atomsDistance,orbitalA,orbitalC);
        }else{
          // SS|SP ; SS|PP
          if (pairTypeB == 1) {
            // SS|SP
            integralValue = IntegralTypeSS_SP(atomsDistance,orbitalA,orbitalC,\
                                              AOsTypeInt,rotationMatrix);
          }else{
            // SS|PP
            integralValue = IntegralTypeSS_PP(atomsDistance,orbitalA,orbitalC,\
                                              AOsTypeInt,rotationMatrix);
          }
        }
      }else{
        // SP|SP ; SP|PP ; PP|PP
        if (pairTypeA == pairTypeB) {
          // SP|SP ; PP|PP
          if (pairTypeA == 1) {
            // SP|SP
            integralValue = IntegralTypeSP_SP(atomsDistance,orbitalA,orbitalC,\
                                              AOsTypeInt,rotationMatrix);
          }else{
            // PP|PP
            integralValue = IntegralTypePP_PP(atomsDistance,orbitalA,orbitalC,\
                                              AOsTypeInt,rotationMatrix);
          }
        }else{
          // SP|PP
            integralValue = IntegralTypeSP_PP(atomsDistance,orbitalA,orbitalC,\
                                              AOsTypeInt,rotationMatrix);
        }
      }
      return convertHartree2eV(integralValue) * invertIntegral;
    }
  }
  return integralValue;
}
/***************************************************************************************/ 
void TwoCenterIntegral::ComputeAllTwoCenterIntegral(const vector<AtomicOrbital>& infoAOs,\
    double**** &all2CenterIntegral){

  int totalAOs = infoAOs.size();
  int indexOfAtomA,indexOfAtomB;
  for (int i=0;i < totalAOs;++i) {
    for (int j=0;j <= i;++j) {
      indexOfAtomA = infoAOs[i].indexAtom;
      indexOfAtomB = infoAOs[j].indexAtom;
      if (infoAOs[i].element == infoAOs[j].element ) {
        // Compute for X-X or H-H
        if (infoAOs[i].element == 1) {
          //Compute for H-H
          ComputePair_HH(i,j,infoAOs,\
              all2CenterIntegral[indexOfAtomA][indexOfAtomB]);
        }else{
          //Compute for X-X
          ComputePair_XX(i,j,infoAOs,\
              all2CenterIntegral[indexOfAtomA][indexOfAtomB]);
        }
      }else{
        // Compute for H-X
        ComputePair_HX(i,j,infoAOs,\
            all2CenterIntegral[indexOfAtomA][indexOfAtomB]);
      }
    }
    if (infoAOs[i].element > 1) {
      i+=3;
    }
  }
}
/***************************************************************************************/ 
void TwoCenterIntegral::ComputePair_HH(const int& indexAO_A,int& indexAO_B,\
    const vector<AtomicOrbital>& infoAOs,double** &all2CenterIntegral){

  all2CenterIntegral[0][0] = ComputeTwoCenterIntegral(infoAOs[indexAO_A],infoAOs[indexAO_A],\
      infoAOs[indexAO_B],infoAOs[indexAO_B]);
}
/***************************************************************************************/ 
void TwoCenterIntegral::ComputePair_HX(int& indexAO_A,int& indexAO_B,\
    const vector<AtomicOrbital>& infoAOs,double** &all2CenterIntegral){
  // List of order lik MOPAC print 
  // SS	SX	XX	SY	XY	YY	SZ	XZ	YZ	ZZ
  int listMOPACA[10] = {0,0,1,0,1,2,0,1,2,3};
  int listMOPACB[10] = {0,1,1,2,2,2,3,3,3,3};
  
  int ajustIndexA,ajustIndexB,ajustIndexC,ajustIndexD;

  bool swapIndexs = false;

  if (infoAOs[indexAO_B].element == 1) {
    swapValues(indexAO_A,indexAO_B);
    swapIndexs = true;
  }

  ajustIndexA = indexAO_A;
  ajustIndexB = indexAO_A;
  
  for (int i=0;i<10;i++) {

    ajustIndexC = listMOPACA[i] + indexAO_B;
    ajustIndexD = listMOPACB[i] + indexAO_B;

    all2CenterIntegral[0][i] = ComputeTwoCenterIntegral(\
        infoAOs[ajustIndexA],infoAOs[ajustIndexB],\
        infoAOs[ajustIndexC],infoAOs[ajustIndexD]);
  }
  
  indexAO_B += 3;

  if (swapIndexs){
    swapValues(indexAO_A,indexAO_B);
  }
}
/***************************************************************************************/ 
void TwoCenterIntegral::ComputePair_XX(const int& indexAO_A,int& indexAO_B,\
       const vector<AtomicOrbital>& infoAOs,double** &all2CenterIntegral){
  // List of order lik MOPAC print 
  // SS	SX	XX	SY	XY	YY	SZ	XZ	YZ	ZZ
  int listMOPACA[10] = {0,0,1,0,1,2,0,1,2,3};
  int listMOPACB[10] = {0,1,1,2,2,2,3,3,3,3};

  int ajustIndexA,ajustIndexB,ajustIndexC,ajustIndexD;

  for (int i=0;i<10;i++) {
    for (int j=0;j<10;j++) {

      ajustIndexA = listMOPACA[i] + indexAO_A;
      ajustIndexB = listMOPACB[i] + indexAO_A;
      ajustIndexC = listMOPACA[j] + indexAO_B;
      ajustIndexD = listMOPACB[j] + indexAO_B;
      all2CenterIntegral[i][j] = ComputeTwoCenterIntegral(\
                                   infoAOs[ajustIndexA],infoAOs[ajustIndexB],\
                                   infoAOs[ajustIndexC],infoAOs[ajustIndexD]);
    }
  }
  indexAO_B += 3;
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
double TwoCenterIntegral::ApplyRotationAOs(int locA, const int& gloA,const double& value,\
    const double (&rotMat)[3][3]){
  return value * rotMat[locA][gloA-1];
}
/***************************************************************************************/ 
double TwoCenterIntegral::ApplyRotationAOs(int locA,int locB,const int& gloA,const int& gloB,\
    const double& value,const double (&rotMat)[3][3]){
  return value * rotMat[locA][gloA-1] * rotMat[locB][gloB-1];
}
/***************************************************************************************/ 
double TwoCenterIntegral::ApplyRotationAOs(int locA,int locB,int locC,\
    const int& gloA,const int& gloB,const int& gloC,\
    const double& value,const double (&rotMat)[3][3]){
  return value * rotMat[locA][gloA-1] * rotMat[locB][gloB-1] * rotMat[locC][gloC-1];
}
/***************************************************************************************/ 
double TwoCenterIntegral::ApplyRotationAOs(int locA,int locB,int locC,int locD,\
    const int& gloA,const int& gloB,const int& gloC, const int& gloD,\
    const double& value,const double (&rotMat)[3][3]){
  return value * rotMat[locA][gloA-1] * rotMat[locB][gloB-1] *\
         rotMat[locC][gloC-1] * rotMat[locD][gloD-1];
}
/***************************************************************************************/ 
double TwoCenterIntegral::IntegralTypeSS_SS(const double& atomsDistance,\
    const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalC){
  return multipole->Interaction_qq(atomsDistance,orbitalA,orbitalC);
}
/***************************************************************************************/ 
double TwoCenterIntegral::IntegralTypeSS_SP(const double& atomsDistance,\
    const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalC,\
    const int (&AOsTypeInt)[4],double (&rotMat)[3][3]){
 
  ComputeRotationMatrix(orbitalA.coordinates,orbitalC.coordinates,rotMat);
  double value = multipole->Interaction_qUz(atomsDistance,orbitalA,orbitalC);
  value = ApplyRotationAOs(2,AOsTypeInt[3],value,rotMat);
  return value;
}
/***************************************************************************************/ 
double TwoCenterIntegral::IntegralTypeSS_PP(const double& atomsDistance,\
    const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalC,\
    const int (&AOsTypeInt)[4],double (&rotMat)[3][3]){
  
  ComputeRotationMatrix(orbitalA.coordinates,orbitalC.coordinates,rotMat);
  
  double value_qq = multipole->Interaction_qq(atomsDistance,orbitalA,orbitalC);
  double value_SS_PpiPpi = value_qq + \
                           multipole->Interaction_qQpipi(atomsDistance,orbitalA,orbitalC);
  value_SS_PpiPpi = ApplyRotationAOs(0,0,AOsTypeInt[2],AOsTypeInt[3],value_SS_PpiPpi,rotMat)\
                    + ApplyRotationAOs(1,1,AOsTypeInt[2],AOsTypeInt[3],value_SS_PpiPpi,rotMat);

  double value_SS_PzPz = value_qq + \
                         multipole->Interaction_qQzz(atomsDistance,orbitalA,orbitalC);
  value_SS_PzPz = ApplyRotationAOs(2,2,AOsTypeInt[2],AOsTypeInt[3],value_SS_PzPz,rotMat);
  return value_SS_PpiPpi + value_SS_PzPz;
}
/***************************************************************************************/ 
double TwoCenterIntegral::IntegralTypeSP_SP(const double& atomsDistance,\
    const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalC,\
    const int (&AOsTypeInt)[4],double (&rotMat)[3][3]){

  ComputeRotationMatrix(orbitalA.coordinates,orbitalC.coordinates,rotMat);
  
  double value_SPpi_SPpi = multipole->Interaction_UpiUpi(atomsDistance,orbitalA,orbitalC);
  ApplyRotationAOs(0,0,AOsTypeInt[1],AOsTypeInt[3],value_SPpi_SPpi,rotMat);

  double value_SPz_SPz = multipole->Interaction_UzUz(atomsDistance,orbitalA,orbitalC);
  ApplyRotationAOs(2,2,AOsTypeInt[1],AOsTypeInt[3],value_SPz_SPz,rotMat);

  return value_SPpi_SPpi + value_SPz_SPz;
}
/***************************************************************************************/ 
double TwoCenterIntegral::IntegralTypeSP_PP(const double& atomsDistance,\
    const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalC,\
    const int (&AOsTypeInt)[4],double (&rotMat)[3][3]){

  ComputeRotationMatrix(orbitalA.coordinates,orbitalC.coordinates,rotMat);
 
  double value_Uzq = multipole->Interaction_Uzq(atomsDistance,orbitalA,orbitalC);

  double value_SPz_PpiPpi = value_Uzq + \
                            multipole->Interaction_UzQpipi(atomsDistance,orbitalA,orbitalC);

  value_SPz_PpiPpi = ApplyRotationAOs(2,0,0,AOsTypeInt[1],AOsTypeInt[2],AOsTypeInt[3],\
      value_SPz_PpiPpi,rotMat) + ApplyRotationAOs(2,1,1,AOsTypeInt[1],AOsTypeInt[2],\
        AOsTypeInt[3],value_SPz_PpiPpi,rotMat);

  double value_SPz_PzPz = value_Uzq + \
                          multipole->Interaction_UzQz(atomsDistance,orbitalA,orbitalC);
  value_SPz_PzPz = ApplyRotationAOs(2,2,2,AOsTypeInt[1],AOsTypeInt[2],AOsTypeInt[3],\
      value_SPz_PzPz,rotMat);

  double value_SPpi_PzPpi = multipole->Interaction_UpiQpiz(atomsDistance,orbitalA,orbitalC);

  value_SPpi_PzPpi = ApplyRotationAOs(0,2,0,AOsTypeInt[1],AOsTypeInt[2],AOsTypeInt[3],\
      value_SPpi_PzPpi,rotMat) + ApplyRotationAOs(1,2,1,AOsTypeInt[1],AOsTypeInt[2],\
        AOsTypeInt[3],value_SPpi_PzPpi,rotMat);

  return value_SPz_PpiPpi + value_SPz_PzPz + value_SPpi_PzPpi;
}
/***************************************************************************************/ 
double TwoCenterIntegral::IntegralTypePP_PP(const double& atomsDistance,\
    const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalC,\
    const int (&AOsTypeInt)[4],double (&rotMat)[3][3]){

  ComputeRotationMatrix(orbitalA.coordinates,orbitalC.coordinates,rotMat);

  double value_qq = multipole->Interaction_qq(atomsDistance,orbitalA,orbitalC);
  double value_qQzz = multipole->Interaction_qQzz(atomsDistance,orbitalA,orbitalC);
  double value_Qzzq = multipole->Interaction_Qzzq(atomsDistance,orbitalA,orbitalC);
  
  double value_PzPz_PzPz = value_qq + value_qQzz + value_Qzzq +\
                           multipole->Interaction_QzzQzz(atomsDistance,orbitalA,orbitalC);

  value_PzPz_PzPz = ApplyRotationAOs(2,2,2,2,AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],\
      AOsTypeInt[3],value_PzPz_PzPz,rotMat);

  double value_qQpipi = multipole->Interaction_qQpipi(atomsDistance,orbitalA,orbitalC);
  double value_Qpipiq = multipole->Interaction_Qpipiq(atomsDistance,orbitalA,orbitalC);

  double value_PpiPpi_PpiPpi = value_qq + value_qQpipi + value_Qpipiq + \
                        multipole->Interaction_QpipiQpipi(atomsDistance,orbitalA,orbitalC);

  double value_PxPy_PxPy = value_PpiPpi_PpiPpi;
  
  value_PpiPpi_PpiPpi = ApplyRotationAOs(0,0,0,0,AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],\
      AOsTypeInt[3],value_PpiPpi_PpiPpi,rotMat) + ApplyRotationAOs(1,1,1,1,\
        AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],AOsTypeInt[3],value_PpiPpi_PpiPpi,rotMat);

  double value_PzPz_PpiPpi = value_qq + value_qQpipi + value_Qzzq + \
                             multipole->Interaction_QzzQpipi(atomsDistance,orbitalA,orbitalC);

  value_PzPz_PpiPpi = ApplyRotationAOs(0,0,0,0,AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],\
      AOsTypeInt[3],value_PzPz_PpiPpi,rotMat) + ApplyRotationAOs(1,1,1,1,\
        AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],AOsTypeInt[3],value_PzPz_PpiPpi,rotMat);

  double value_PpiPpi_PzPz = value_qq + value_Qpipiq + value_qQzz + \
                             multipole->Interaction_QpipiQzz(atomsDistance,orbitalA,orbitalC);

  value_PpiPpi_PzPz = ApplyRotationAOs(0,0,0,0,AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],\
      AOsTypeInt[3],value_PpiPpi_PzPz,rotMat) + ApplyRotationAOs(1,1,1,1,\
        AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],AOsTypeInt[3],value_PpiPpi_PzPz,rotMat);

  double value_PpiPz_PpiPz = multipole->Interaction_QpizQpiz(atomsDistance,orbitalA,orbitalC); 

  value_PpiPz_PpiPz = ApplyRotationAOs(0,0,0,0,AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],\
      AOsTypeInt[3],value_PpiPz_PpiPz,rotMat) + ApplyRotationAOs(1,1,1,1,\
        AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],AOsTypeInt[3],value_PpiPz_PpiPz,rotMat);

  double value_PxPx_PyPy = value_qq + value_qQpipi + value_Qpipiq + \
                        multipole->Interaction_QxxQyy(atomsDistance,orbitalA,orbitalC);
  
  value_PxPy_PxPy = (value_PxPy_PxPy - value_PxPx_PyPy) * 0.5;

  value_PxPx_PyPy = ApplyRotationAOs(0,0,1,1,AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],\
      AOsTypeInt[3],value_PxPx_PyPy,rotMat) + ApplyRotationAOs(1,1,0,0,\
        AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],AOsTypeInt[3],value_PxPx_PyPy,rotMat);

  value_PxPy_PxPy = ApplyRotationAOs(0,1,0,1,AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],\
      AOsTypeInt[3],value_PxPx_PyPy,rotMat) + ApplyRotationAOs(1,1,0,0,\
        AOsTypeInt[0],AOsTypeInt[1],AOsTypeInt[2],AOsTypeInt[3],value_PxPx_PyPy,rotMat);

  return value_PzPz_PzPz + value_PpiPpi_PpiPpi + value_PzPz_PpiPpi + value_PpiPpi_PzPz + \
    value_PpiPz_PpiPz + value_PxPx_PyPy + value_PxPy_PxPy;
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
double TwoCenterIntegral::CorrectOrderIntegral(int& pairTypeA,int& pairTypeB,\
    int (&AOsTypeInt)[4]){
  if (pairTypeA > pairTypeB) {
    int tmp_pair = pairTypeA;
    pairTypeA = pairTypeB;
    pairTypeB = tmp_pair;

    int AOsTypeInt_tmp = AOsTypeInt[0];
    int AOsTypeInt_tmpp = AOsTypeInt[1];

    AOsTypeInt[0] = AOsTypeInt[2];
    AOsTypeInt[1] = AOsTypeInt[3];

    AOsTypeInt[2] = AOsTypeInt_tmp;
    AOsTypeInt[3] = AOsTypeInt_tmpp;

    if ((pairTypeA + pairTypeB) % 2 == 1) {
      return -1.0;
    }else{
      return 1.0;
    }
  }
  return 1.0;
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _TWOCENTERINTEGRAL_H_
