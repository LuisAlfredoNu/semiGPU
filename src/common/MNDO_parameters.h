/* 
   List of parameters for semiempiric method MNDO
*/
#ifndef _MNDO_PARAMETERS_H_
#define _MNDO_PARAMETERS_H_

/***************************************************************************************/ 
class MNDOparameter{
/***************************************************************************************/ 
public:
   MNDOparameter();
   double uss[109];
   double upp[109];
   
   double zs[109];

   double betaS[109];
   double betaP[109];

   double alpha[109];

   double gss[109];
   double gpp[109];
   double gsp[109];
   double gp2[109];
   double hsp[109];
//TODO Is posible compute this values 
//  ref: https://arxiv.org/pdf/1806.06147v3
   double dd[109][2];
   double pp[109][3];
//TODO END

   double heat[109];
   double eisol[109];
};

#endif // _MNDO_PARAMETERS_H_
