/* 
   List of parameters for semiempiric method MNDO
*/
#ifndef _MNDO_PARAMETERS_CPP_
#define _MNDO_PARAMETERS_CPP_

#include "MNDO_parameters.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

MNDOparameter::MNDOparameter(){

// Data for Hidrogen
 
  uss[1] = -11.90627600;  // EV
  zs[1] = 1.3319670;      // a.u.
  betaS[1] = -6.98906400; // EV
  alpha[1] = 2.54413410;  // 1/Ang

  gss[1] = 12.84800000;   // EV
  
  p01[1] = 1.05897351;    // Bohr

  heat[1] = 52.102;       // kcal/mol
  eisol[1] = -11.9062760; // EV

// Data for Carbon
  uss[6] = -52.27974500;  // EV
  upp[6] = -39.20555800;  // EV
  zs[6] = 1.7875370;      // a.u.
  betaS[6] = -18.9850440; // EV
  betaP[6] = -7.93412200; // EV
  alpha[6] = 2.546380;    // 1/Ang

  gss[6] = 12.23000000;   // EV
  gpp[6] = 11.08000000;   // EV
  gsp[6] = 11.47000000;   // EV
  gp2[6] = 9.84000000;    // EV
  hsp[6] = 2.43000000;    // EV

  dd2[6] = 0.80746618;    // Bohr
  dd3[6] = 0.68515777;    // Bohr
  p01[6] = 1.11248501;    // Bohr
  p02[6] = 0.81309641;    // Bohr
  p03[6] = 0.74785493;    // Bohr

  heat[6] = 170.89;       // kcal/mol
  eisol[6] = -120.500606; // EV

}



#endif // _MNDO_PARAMETERS_CPP_
