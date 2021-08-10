/*
 * Class for compute the multipole interactions   
 */

#ifndef _MULTIPOLE_CPP_
#define _MULTIPOLE_CPP_

#include "atomicOrbitals.h"

#include <cmath>

#include "multipole.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

Multipole::Multipole(MNDOparameter& MNDOpara){
  parameter = &MNDOpara ;
}
/***************************************************************************************/ 
double Multipole::Interaction_qq(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){

  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,0,orbitalB.element,0);
  double integral;
  integral = sqrt(atomDistance * atomDistance + addTerm);
  integral = 1.0 / integral;
  
  return integral;
}
/***************************************************************************************/ 
double Multipole::Interaction_qUz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,0,orbitalB.element,1);
  double finalIntegral = 0.0e-12;
  double integral;
  integral = atomDistance + parameter->dd[orbitalB.element][0] ;
  integral =  2.0 * sqrt(integral * integral + addTerm);
  integral = 1.0 / integral;
  
  finalIntegral += integral; 
  
  integral = atomDistance - parameter->dd[orbitalB.element][0] ;
  integral =  2.0 * sqrt(integral * integral + addTerm);
  integral = 1.0 / integral;
  
  finalIntegral -= integral;
  
  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_Uzq(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  return Interaction_qUz(atomDistance,orbitalB,orbitalA) * -1.0 ;
}
/***************************************************************************************/ 
double Multipole::Interaction_qQpipi(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,0,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;

  integral = std::pow(2.0 * parameter->dd[orbitalB.element][1],2.0);
  integral = 2.0 * sqrt(atomDistance * atomDistance + integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral += integral;

  integral = 2.0 * sqrt(atomDistance * atomDistance + addTerm);
  integral = 1.0 / integral;

  finalIntegral -= integral;

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_Qpipiq(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  return Interaction_qQpipi(atomDistance,orbitalB,orbitalA) ;
}
/***************************************************************************************/ 
double Multipole::Interaction_qQzz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,0,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;

  integral = atomDistance + 2.0 * parameter->dd[orbitalB.element][1];
  integral = 4.0 * sqrt(integral * integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral += integral;

  integral = 2.0 * sqrt(atomDistance * atomDistance + addTerm);
  integral = 1.0 / integral;

  finalIntegral -= integral;
  
  integral = atomDistance - 2.0 * parameter->dd[orbitalB.element][1];
  integral = 4.0 * sqrt(integral * integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral += integral;

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_Qzzq(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  return Interaction_qQzz(atomDistance,orbitalB,orbitalA) ;
}
/***************************************************************************************/ 
double Multipole::Interaction_UpiUpi(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,1,orbitalB.element,1);
  double finalIntegral = 0.0e-12;
  double integral;
  
  integral = parameter->dd[orbitalA.element][0] - parameter->dd[orbitalB.element][0];
  integral = 2.0 * sqrt(atomDistance * atomDistance + integral * integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral += integral;
  
  integral = parameter->dd[orbitalA.element][0] + parameter->dd[orbitalB.element][0];
  integral = 2.0 * sqrt(atomDistance * atomDistance + integral * integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral -= integral;
  
  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_UzUz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,1,orbitalB.element,1);
  double finalIntegral = 0.0e-12;
  double integral;
  double swapSign[2] = {1.0,-1.0};
  double swapSignn[3] = {-1.0,1.0,-1.0};

  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      integral = swapSign[i] * parameter->dd[orbitalA.element][0] \
                 +  swapSign[j] * parameter->dd[orbitalB.element][0] + atomDistance;
      integral = 4.0 * sqrt(integral * integral + addTerm);
      integral = 1.0 / integral;
      
      finalIntegral += swapSignn[i+j] * integral;
    }
  }

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_UpiQpiz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,1,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;
  double swapSign[2] = {1.0,-1.0};
  double swapSignn[3] = {-1.0,1.0,-1.0};

  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      integral = std::pow(atomDistance \
                 + swapSign[i] * parameter->dd[orbitalB.element][1],2.0) \
                 + std::pow(parameter->dd[orbitalA.element][0] \
                    + swapSign[j] * parameter->dd[orbitalB.element][1],2.0);
      integral = 4.0 * sqrt(integral + addTerm);
      integral = 1.0 / integral;
      
      finalIntegral += swapSignn[i+j] * integral;
    }
  }

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_QpizUpi(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  return Interaction_UpiQpiz(atomDistance,orbitalB,orbitalA) * -1.0;
}
/***************************************************************************************/ 
double Multipole::Interaction_UzQpipi(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,1,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;
  double swapSign[2] = {1.0,-1.0};
  double swapSignn[2] = {-1.0,1.0};

  for (int i=0;i<2;i++) {
    integral = std::pow(atomDistance \
        + swapSign[i] * parameter->dd[orbitalA.element][0],2.0) \
               + std::pow(2.0 * parameter->dd[orbitalB.element][1],2.0);
    integral = 4.0 * sqrt(integral + addTerm);
    integral = 1.0 / integral;

    finalIntegral += swapSignn[i] * integral;
  }
  for (int i=0;i<2;i++) {
    integral = std::pow(atomDistance \
        + swapSign[i] * parameter->dd[orbitalA.element][0],2.0);
    integral = 4.0 * sqrt(integral + addTerm);
    integral = 1.0 / integral;

    finalIntegral += swapSign[i] * integral;
  }

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_QpipiUz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  return Interaction_UzQpipi(atomDistance,orbitalB,orbitalA) * -1.0;
}
/***************************************************************************************/ 
double Multipole::Interaction_UzQz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,1,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;
  double swapSign[2] = {1.0,-1.0};
  double swapSignn[4] = {-1.0,-1.0,1.0,1.0};

  int count=0;

  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      integral = atomDistance \
          + swapSign[i] * parameter->dd[orbitalA.element][0] \
          + swapSign[j] * 2.0 * parameter->dd[orbitalB.element][1];
      integral = 8.0 * sqrt(integral * integral + addTerm);
      integral = 1.0 / integral;

      finalIntegral += swapSignn[count++] * integral;
    }
  }
  for (int i=0;i<2;i++) {
    integral = atomDistance + swapSign[i] * parameter->dd[orbitalA.element][0];
    integral = 4.0 * sqrt(integral * integral + addTerm);
    integral = 1.0 / integral;

    finalIntegral += swapSign[i] * integral;
  }

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_QzUz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  return Interaction_UzQz(atomDistance,orbitalB,orbitalA) * -1.0;
}
/***************************************************************************************/ 
double Multipole::Interaction_QpipiQpipi(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,2,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;
  double swapSign[2] = {1.0,-1.0};

  for (int i=0;i<2;i++) {
    integral = 4.0 * std::pow( parameter->dd[orbitalA.element][1] \
               + swapSign[i] * parameter->dd[orbitalB.element][1],2.0);
    integral = 8.0 * sqrt(atomDistance * atomDistance + integral + addTerm);
    integral = 1.0 / integral;

    finalIntegral +=  integral;
  }
  integral = std::pow(2.0 * parameter->dd[orbitalA.element][1],2.0);
  integral = 4.0 * sqrt(atomDistance * atomDistance + integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral -=  integral;

  integral = std::pow(2.0 * parameter->dd[orbitalB.element][1],2.0);
  integral = 4.0 * sqrt(atomDistance * atomDistance + integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral -=  integral;

  integral = 4.0 * sqrt(atomDistance * atomDistance + addTerm);
  integral = 1.0 / integral;

  finalIntegral +=  integral;

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_QxxQyy(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,2,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;

  integral =  std::pow( 2.0 * parameter->dd[orbitalA.element][1],2.0) \
              + std::pow( 2.0 * parameter->dd[orbitalB.element][1],2.0);
  integral = 4.0 * sqrt(atomDistance * atomDistance + integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral +=  integral;

  integral = std::pow(2.0 * parameter->dd[orbitalA.element][1],2.0);
  integral = 4.0 * sqrt(atomDistance * atomDistance + integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral -=  integral;

  integral = std::pow(2.0 * parameter->dd[orbitalB.element][1],2.0);
  integral = 4.0 * sqrt(atomDistance * atomDistance + integral + addTerm);
  integral = 1.0 / integral;

  finalIntegral -=  integral;

  integral = 4.0 * sqrt(atomDistance * atomDistance + addTerm);
  integral = 1.0 / integral;

  finalIntegral +=  integral;

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_QpipiQzz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,2,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral,integralAux,integralAuxx;
  double swapSign[2] = {1.0,-1.0};

  integralAuxx = std::pow(2.0 * parameter->dd[orbitalA.element][1],2.0) + addTerm;

  for (int i=0;i<2;i++) {

    integralAux = std::pow(atomDistance \
                  + 2.0 *  swapSign[i] * parameter->dd[orbitalB.element][1] ,2.0);

    integral = 8.0 * sqrt(integralAux +integralAuxx);
    integral = 1.0 / integral;

    finalIntegral +=  integral;
    integral = 8.0 * sqrt(integralAux + addTerm);
    integral = 1.0 / integral;

    finalIntegral -=  integral;

    integral = 8.0 * sqrt(atomDistance * atomDistance + integralAuxx);
    integral = 1.0 / integral;

    finalIntegral -=  integral;
  }
  integral = 4.0 * sqrt(atomDistance * atomDistance + addTerm);
  integral = 1.0 / integral;

  finalIntegral +=  integral;
  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_QzzQpipi(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  return Interaction_QpipiQzz(atomDistance,orbitalB,orbitalA);
}
/***************************************************************************************/ 
double Multipole::Interaction_QzzQzz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,2,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;
  double swapSign[2] = {1.0,-1.0};

  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      integral = atomDistance \
                 + swapSign[i] * 2.0 * parameter->dd[orbitalA.element][1] \
                 + swapSign[j] * 2.0 * parameter->dd[orbitalB.element][1];
      integral = 16.0 * sqrt(integral * integral + addTerm);
      integral = 1.0 / integral;

      finalIntegral +=  integral;
    }
  }
  for (int i=0;i<2;i++) {
    integral = atomDistance + 2.0 * swapSign[i] * parameter->dd[orbitalA.element][1];
    integral = 8.0 * sqrt(integral * integral + addTerm);
    integral = 1.0 / integral;

    finalIntegral -= integral;
    
    integral = atomDistance + 2.0 * swapSign[i] * parameter->dd[orbitalB.element][1];
    integral = 8.0 * sqrt(integral * integral + addTerm);
    integral = 1.0 / integral;

    finalIntegral -= integral;
  }
  
  integral = 4.0 * sqrt(atomDistance * atomDistance + addTerm);
  integral = 1.0 / integral;

  finalIntegral +=  integral;

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_QpizQpiz(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,2,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;
  double swapSign[2] = {1.0,-1.0};
  double swapSignn[8] = {1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0};

  int count = 0;

  for (int i=0;i<2;i++) {
    for (int j=0;j<2;j++) {
      for (int k=0;k<2;k++) {
        integral = std::pow( atomDistance \
            + swapSign[i] * parameter->dd[orbitalA.element][1] \
            + swapSign[j] * parameter->dd[orbitalB.element][1],2.0);
        integral = integral + std::pow( parameter->dd[orbitalA.element][1] \
            + swapSign[k] * parameter->dd[orbitalB.element][1],2.0);
        integral = 8.0 * sqrt(integral + addTerm);
        integral = 1.0 / integral;

        finalIntegral +=  swapSignn[count++] * integral;
      }
    }
  }

  return finalIntegral;
}
/***************************************************************************************/ 
double Multipole::Interaction_QxyQxy(const double& atomDistance,const AtomicOrbital& orbitalA,\
    const AtomicOrbital& orbitalB){
  
  double addTerm;
  AditiveTerms(addTerm,orbitalA.element,2,orbitalB.element,2);
  double finalIntegral = 0.0e-12;
  double integral;
  double swapSign[2] = {1.0,-1.0};

  for (int i=0;i<2;i++) {

    integral = parameter->dd[orbitalA.element][1] \
                  + swapSign[i] * parameter->dd[orbitalB.element][1];
    integral *= integral;

    integral = 4.0 * sqrt(atomDistance * atomDistance + 2.0 * integral + addTerm);
    integral = 1.0 / integral;

    finalIntegral +=  integral;

    integral = 2.0 * pow(parameter->dd[orbitalA.element][1],2.0) \
              + 2.0 * pow(parameter->dd[orbitalB.element][1],2.0);
    
    integral = 4.0 * sqrt(atomDistance * atomDistance + integral + addTerm);
    integral = 1.0 / integral;

    finalIntegral -=  integral;
  }

  return finalIntegral;
}
/***************************************************************************************/ 
void Multipole::AditiveTerms(double& addTerm,const int& elementA,int l,\
       const int& elementB,int m){
  addTerm = parameter->pp[elementA][l] + parameter->pp[elementB][m];
  addTerm *= addTerm;
}

#endif // _MULTIPOLE_CPP_
