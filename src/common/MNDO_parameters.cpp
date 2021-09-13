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
  
  pp[1][0] = 1.05897351;    // Bohr

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

  dd[6][0] = 0.80746618;    // Bohr
  dd[6][1] = 0.68515777;    // Bohr
  pp[6][0] = 1.11248501;    // Bohr
  pp[6][1] = 0.81309641;    // Bohr
  pp[6][2] = 0.74785493;    // Bohr

   heat[6] = 170.89;       // kcal/mol
  eisol[6] = -120.500606; // EV
// Data for Nitrogen
   uss[7] =  -71.93212200; //  EV        ONE-CENTER ENERGY FOR S
   upp[7] =  -57.17231900; //  EV        ONE-CENTER ENERGY FOR P
    zs[7] =    2.25561400; //  AU        ORBITAL EXPONENT  FOR S
 betaS[7] =  -20.49575800; //  EV        BETA PARAMETER    FOR S
 betaP[7] =  -20.49575800; //  EV        BETA PARAMETER    FOR P
 alpha[7] =    2.86134200; //  (1/A)     ALPHA PARAMETER   FOR CORE
   gss[7] =   13.59000000; //  EV        ONE-CENTER INTEGRAL (SS,SS)
   gpp[7] =   12.98000000; //  EV        ONE-CENTER INTEGRAL (PP,PP)
   gsp[7] =   12.66000000; //  EV        ONE-CENTER INTEGRAL (SS,PP)
   gp2[7] =   11.59000000; //  EV        ONE-CENTER INTEGRAL (PP*,PP*)
   hsp[7] =    3.14000000; //  EV        ONE-CENTER INTEGRAL (SP,SP)
 dd[7][0] =    0.63990367; //  BOHR      CHARGE SEPARATION, SP, L=1
 dd[7][1] =    0.76788440; //  BOHR      CHARGE SEPARATION, PP, L=2
 dd[7][1] =    0.54297627; //  BOHR      USING ORIGINAL MNDO PAPER FORMULA
 pp[7][0] =    1.00115465; //  BOHR      KLOPMAN-OHNO TERM, SS, L=0
 pp[7][1] =    0.63747362; //  BOHR      KLOPMAN-OHNO TERM, SP, L=1
 pp[7][2] =    0.61528499; //  BOHR      KLOPMAN-OHNO TERM, PP, L=2
  heat[7] =  113.00000000; //  KCAL/MOL  HEAT OF FORMATION OF THE ATOM (EXP)
 eisol[7] = -202.56620100; //  EV        TOTAL ENERGY OF THE ATOM (CALC)

// Data for Oxygen 
   uss[8] =  -99.64430900  ; // EV        ONE-CENTER ENERGY FOR S
   upp[8] =  -77.79747200  ; // EV        ONE-CENTER ENERGY FOR P
    zs[8] =    2.69990500  ; // AU        ORBITAL EXPONENT  FOR S
 betaS[8] =   -32.68808200 ; //  EV        BETA PARAMETER    FOR S
 betaP[8] =   -32.68808200 ; //  EV        BETA PARAMETER    FOR P
 alpha[8] =     3.16060400 ; //  (1/A)     ALPHA PARAMETER   FOR CORE
   gss[8] =    15.42000000 ; //  EV        ONE-CENTER INTEGRAL (SS,SS)
   gpp[8] =    14.52000000 ; //  EV        ONE-CENTER INTEGRAL (PP,PP)
   gsp[8] =    14.48000000 ; //  EV        ONE-CENTER INTEGRAL (SS,PP)
   gp2[8] =    12.98000000 ; //  EV        ONE-CENTER INTEGRAL (PP*,PP*)
   hsp[8] =     3.94000000 ; //  EV        ONE-CENTER INTEGRAL (SP,SP)
 dd[8][0] =     0.53460239 ; //  BOHR      CHARGE SEPARATION, SP, L=1
 dd[8][1] =     0.64152287 ; //  BOHR      CHARGE SEPARATION, PP, L=2
 dd[8][1] =     0.45362517 ; //  BOHR      USING ORIGINAL MNDO PAPER FORMULA
 pp[8][0] =     0.88234058 ; //  BOHR      KLOPMAN-OHNO TERM, SS, L=0
 pp[8][1] =     0.52124931 ; //  BOHR      KLOPMAN-OHNO TERM, SP, L=1
 pp[8][2] =     0.52654944 ; //  BOHR      KLOPMAN-OHNO TERM, PP, L=2
  heat[8] =    59.55900000 ; //  KCAL/MOL  HEAT OF FORMATION OF THE ATOM (EXP)
 eisol[8] =  -317.86850600 ; //  EV        TOTAL ENERGY OF THE ATOM (CALC)


}



#endif // _MNDO_PARAMETERS_CPP_
