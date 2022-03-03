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
   ~MNDOparameter();

/***************************************************************************************/
// Variables
/***************************************************************************************/ 
   unsigned int countParameter_;

   double* uss;
   double* upp;
   
   double* zs;

   double* betaS;
   double* betaP;

   double* alpha;

   double* gss;
   double* gpp;
   double* gsp;
   double* gp2;
   double* hsp;
//TODO Is posible compute this values 
//  ref: https://arxiv.org/pdf/1806.06147v3
   double** dd;//[2]
   double** pp;//[3]
//TODO END

   double* heat;
   double* eisol;
/***************************************************************************************/ 
/***************************************************************************************/ 
private:
/***************************************************************************************/ 
// Methods
/***************************************************************************************/ 
   void SetArraysSizes();
   void SetParameterValues();
   void To_device();
   void From_device();
};

#endif // _MNDO_PARAMETERS_H_
