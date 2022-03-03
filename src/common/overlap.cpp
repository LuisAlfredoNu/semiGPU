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
#include "STO-6G.h"

#include "overlap.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
Overlap::Overlap(const STO_6G& basisSTO_) {
  basisSTO = &basisSTO_;
  xk_ = new double[12];
  wk_ = new double[12];

  xk_[0]  = -3.889724897869781919272;       xk_[1]  = -3.020637025120889771711;
  xk_[2]  = -2.279507080501059900188;       xk_[3]  = -1.59768263515260479671;
  xk_[4]  = -0.9477883912401637437046;      xk_[5]  = -0.314240376254359111277;
  xk_[6]  =  3.889724897869781919272;       xk_[7]  =  3.020637025120889771711;
  xk_[8]  =  2.279507080501059900188;       xk_[9]  =  1.59768263515260479671;
  xk_[10] =  0.9477883912401637437046;      xk_[11] =  0.314240376254359111277;
  wk_[0]  =  2.65855168435630160602E-7;     wk_[1]  =  8.5736870435878586546E-5;
  wk_[2]  =  0.00390539058462906185999;     wk_[3]  =  0.05160798561588392999187;
  wk_[4]  =  0.2604923102641611292334;      wk_[5]  =  0.5701352362624795783471;
  wk_[6]  =  2.65855168435630160602E-7;     wk_[7]  =  8.5736870435878586546E-5;
  wk_[8]  =  0.00390539058462906185999;     wk_[9]  =  0.05160798561588392999187;
  wk_[10] =  0.2604923102641611292334;      wk_[11] =  0.5701352362624795783471;
  #pragma acc enter data copyin(this[0:1],this->xk_[0:12],this->wk_[0:12])
}
/***************************************************************************************/ 
Overlap::Overlap(const ListAtomicOrbitals &infoAOs,const STO_6G& basisSTO_) : BaseMatrix(infoAOs.size()){
  basisSTO = &basisSTO_;
  infoAOs_ = &infoAOs;
  xk_ = new double[12];
  wk_ = new double[12];

  xk_[0]  = -3.889724897869781919272;       xk_[1]  = -3.020637025120889771711;
  xk_[2]  = -2.279507080501059900188;       xk_[3]  = -1.59768263515260479671;
  xk_[4]  = -0.9477883912401637437046;      xk_[5]  = -0.314240376254359111277;
  xk_[6]  =  3.889724897869781919272;       xk_[7]  =  3.020637025120889771711;
  xk_[8]  =  2.279507080501059900188;       xk_[9]  =  1.59768263515260479671;
  xk_[10] =  0.9477883912401637437046;      xk_[11] =  0.314240376254359111277;
  wk_[0]  =  2.65855168435630160602E-7;     wk_[1]  =  8.5736870435878586546E-5;
  wk_[2]  =  0.00390539058462906185999;     wk_[3]  =  0.05160798561588392999187;
  wk_[4]  =  0.2604923102641611292334;      wk_[5]  =  0.5701352362624795783471;
  wk_[6]  =  2.65855168435630160602E-7;     wk_[7]  =  8.5736870435878586546E-5;
  wk_[8]  =  0.00390539058462906185999;     wk_[9]  =  0.05160798561588392999187;
  wk_[10] =  0.2604923102641611292334;      wk_[11] =  0.5701352362624795783471;
}
/***************************************************************************************/ 
double Overlap::ComputeElementMatrix(const size_t &i,const size_t &j){
  //return ComputeOverlap_Boys(infoAOs_->orbital[(int)i],infoAOs_->orbital[(int)j]);
  return ComputeOverlap(infoAOs_->orbital[i],infoAOs_->orbital[j]);
}
/***************************************************************************************/ 
double Overlap::ComputeOverlap(const AtomicOrbital& orbitalA,\
                               const AtomicOrbital& orbitalB){
  // Code develop by Julio Céar Cruz Monterrosas Feb 2021.

  double AB = 0.0e-10;
  for (int i=0;i<3;++i) {
    AB += (orbitalA.coordinates[i] - orbitalB.coordinates[i]) * \
          (orbitalA.coordinates[i] - orbitalB.coordinates[i]);
  }
  int sumAngularMomentumA = orbitalA.angularMomentum[0] + orbitalA.angularMomentum[1] + \
                            orbitalA.angularMomentum[2];
  int sumAngularMomentumB = orbitalB.angularMomentum[0] + orbitalB.angularMomentum[1] + \
                            orbitalB.angularMomentum[2];
  double P,PA[3],PB[3],mu, p;
  double alpha, beta;
  double overlapTot,overlapFinal = 0.0e-10;
  for (int gtoA=0;gtoA<6;++gtoA) {
    for (int gtoB=0;gtoB<6;++gtoB) {

      if (sumAngularMomentumA == 0) {
        alpha = basisSTO->value[orbitalA.element].exponentS[gtoA];
      }else{
        alpha = basisSTO->value[orbitalA.element].exponentP[gtoA];
      }
      if (sumAngularMomentumB == 0) {
        beta  = basisSTO->value[orbitalB.element].exponentS[gtoB];
      }else{
        beta  = basisSTO->value[orbitalB.element].exponentP[gtoB];
      }
      p = alpha + beta;
      mu = alpha*beta/p;
  
      overlapTot = exp(-mu * AB);
     
      if (sumAngularMomentumA == 0) {
        overlapTot *= basisSTO->value[orbitalA.element].coeffS[gtoA];
      }else{
        overlapTot *= basisSTO->value[orbitalA.element].coeffP[gtoA];
      }
      if (sumAngularMomentumB == 0) {
        overlapTot *= basisSTO->value[orbitalB.element].coeffS[gtoB];
      }else{
        overlapTot *= basisSTO->value[orbitalB.element].coeffP[gtoB];
      }
      //cout << "preExp = " << overlapTot << endl;

      for (int i=0; i<3; ++i) {
        P  = alpha*orbitalA.coordinates[i];
        P +=  beta*orbitalB.coordinates[i];
        P /= p;
        PA[i] = P - orbitalA.coordinates[i];
        PB[i] = P - orbitalB.coordinates[i];
        //overlapTot *= OverlapMcMurchie(orbitalA.angularMomentum[i],\
                       orbitalB.angularMomentum[i], p, PA[i], PB[i]);
        overlapTot *= OverlapNumericalIntegral(orbitalA.angularMomentum[i],\
                       orbitalB.angularMomentum[i], p, PA[i], PB[i]);
      }
     /** 
      overlapTot *= NormalizationConst(alpha,orbitalA.angularMomentum);
      overlapTot *= NormalizationConst(beta,orbitalB.angularMomentum);
      
      **/
      overlapFinal += overlapTot ;
    }
  }
  return overlapFinal;
}
/***************************************************************************************/ 
double Overlap::McMurchieCoef (const int &la,const int &lb,const int &t,const double &p,const double &Xpa,const double &Xpb){
  if (t < 0 || t > la+lb)
    return 0.0;
  if (t == 0) {
    if (la == 0 && lb == 0) {
      return 1.0;
    }else if (la > lb){
        return Xpa * McMurchieCoef(la-1, lb, 0, p, Xpa, Xpb) +\
                     McMurchieCoef(la-1, lb, 1, p, Xpa, Xpb);
    }else{
      return Xpb * McMurchieCoef(la, lb-1, 0, p, Xpa, Xpb) + \
                   McMurchieCoef(la, lb-1, 1, p, Xpa, Xpb);
    }
  } else if (t == 1){
    return (1.0 / (2.0 * p)) * (la * McMurchieCoef(la-1, lb, 0, p, Xpa, Xpb) + \
                                lb * McMurchieCoef(la, lb-1, 0, p, Xpa, Xpb));
  }else{
    return (1.0 / (2.0 * p * t)) * (la * McMurchieCoef(la-1, lb, t-1, p, Xpa, Xpb) +\
                                    lb * McMurchieCoef(la, lb-1, t-1, p, Xpa, Xpb));
  }
}
/***************************************************************************************/ 
double Overlap::OverlapMcMurchie(const int &la,const int &lb,const double &expoTot,const double &Xpa,const double &Xpb) {
  if (la == 0 && lb == 0){
    return std::sqrt(M_PI/expoTot);
  }else{
    return std::sqrt(M_PI/expoTot) * McMurchieCoef(la, lb, 0, expoTot, Xpa, Xpb);
  }
}
/***************************************************************************************/
double Overlap::NormalizationConst(const double &alpha,const int (&angMom)[3]){
  double norm = std::pow(2.0 * alpha/M_PI,3.0/4.0);
  norm *= std::pow(4.0 * alpha,(double)(angMom[0]+angMom[1]+angMom[2])/2);
  norm /= std::sqrt(\
                    doubfact(2*angMom[0]-1) * \
                    doubfact(2*angMom[1]-1) * \
                    doubfact(2*angMom[2]-1)   \
                   );
  return norm;
}
/***************************************************************************************/
double Overlap::OverlapNumericalIntegral(const int &la,const int &lb,const double &expoTot,const double &Xpa,const double &Xpb) {
  
  double sqrt_p = std::sqrt(expoTot);

  if (la == 0 && lb == 0){
    return std::sqrt(M_PI) / sqrt_p;
  }

  double x,integral = 0.0;
  
  for (short i=0;i<12;++i) {
    x = xk_[i]/sqrt_p;
    integral += wk_[i] * pow(x + Xpa, la) * pow(x + Xpb, lb);
  }
  return integral/sqrt_p;
}
/***************************************************************************************/ 
double Overlap::ComputeOverlap_Boys(const AtomicOrbital& orbitalA,const AtomicOrbital& orbitalB){
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
  double alpha,beta;
  double a_times_b = 0.0;
  double a_plus_b = 0.0;
  // Normalization factor
  double normaFactor;

  for(int gtoA = 0; gtoA < 6;gtoA++){
    for(int gtoB = 0; gtoB < 6;gtoB++){
      if (sumAngularMomentumA == 0) {
        alpha = basisSTO->value[orbitalA.element].exponentS[gtoA];
      }else{
        alpha = basisSTO->value[orbitalA.element].exponentP[gtoA];
      }
      if (sumAngularMomentumB == 0) {
        beta  = basisSTO->value[orbitalB.element].exponentS[gtoB];
      }else{
        beta  = basisSTO->value[orbitalB.element].exponentP[gtoB];
      }
      a_times_b = alpha * beta;
      a_plus_b = alpha + beta;

      Overlap_SS(overlapGTO,a_times_b,a_plus_b,atomDistance);

      if (sumAngularMomentumA == 1 && sumAngularMomentumB == 0) {
        Overlap_PS(overlapGTO,alpha,\
                              beta,\
                              a_plus_b,rABx);
      
        normaFactor = basisSTO->value[atomTypeA].coeffP[gtoA] * \
                      basisSTO->value[atomTypeB].coeffS[gtoB];
      }else if (sumAngularMomentumA == 0 && sumAngularMomentumB == 1) {
        Overlap_SP(overlapGTO,alpha,\
                              beta,\
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
