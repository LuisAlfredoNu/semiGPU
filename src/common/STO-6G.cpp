#ifndef _STO_6G_CPP_
#define _STO_6G_CPP_

#include "STO-6G.h"

STO_6G::STO_6G(){
  countAtomValues_ = 10;
  SetArraysSizes();
  SetParameterValues();
  To_device();
}
/***************************************************************************************/ 
STO_6G::~STO_6G(){
  From_device();
}
/***************************************************************************************/ 
void STO_6G::SetArraysSizes(){
  value = new STOdata[countAtomValues_];
}
/***************************************************************************************/ 
void STO_6G::SetParameterValues(){
// Hidrogen 1s
  value[1].exponentS[0] = 0.3552322122e+02;       
  value[1].exponentS[1] = 0.6513143725e+01;       
  value[1].exponentS[2] = 0.1822142904e+01;       
  value[1].exponentS[3] = 0.6259552659e+00;       
  value[1].exponentS[4] = 0.2430767471e+00;       
  value[1].exponentS[5] = 0.1001124280e+00;      
  
     value[1].coeffS[0] = 0.9163596281e-02;
     value[1].coeffS[1] = 0.4936149294e-01;
     value[1].coeffS[2] = 0.1685383049e+00;
     value[1].coeffS[3] = 0.3705627997e+00;
     value[1].coeffS[4] = 0.4164915298e+00;
     value[1].coeffS[5] = 0.1303340841e+00;
  
  /***************************************************************************************/ 
  //
  // Stolen from Sparrow
  //
  /***************************************************************************************/ 
   value[1].exponentS[0] =   4.098792193376860e+01;   
   value[1].exponentS[1] =   7.515090619202102e+00;   
   value[1].exponentS[2] =   2.102451537968084e+00;   
   value[1].exponentS[3] =   7.222488471249888e-01;   
   value[1].exponentS[4] =   2.804703624957924e-01;   
   value[1].exponentS[5] =   1.155131838418812e-01;   
   
      value[1].coeffS[0] =   1.057956874127020e-01;
      value[1].coeffS[1] =   1.596794562489045e-01;
      value[1].coeffS[2] =   2.097264517638662e-01;
      value[1].coeffS[3] =   2.069126553471125e-01;
      value[1].coeffS[4] =   1.144014417359658e-01;
      value[1].coeffS[5] =   1.840525175585220e-02;
  
  /***************************************************************************************/ 
  //
  // End Stolen from Sparrow
  //
  /***************************************************************************************/ 
/***************************************************************************************/ 
// Carbon 2s, 2p
  value[6].exponentS[0] = 0.3049723950e+02;
  value[6].exponentS[1] = 0.6036199601e+01;
  value[6].exponentS[2] = 0.1876046337e+01;
  value[6].exponentS[3] = 0.7217826470e+00;
  value[6].exponentS[4] = 0.3134706954e+00;
  value[6].exponentS[5] = 0.1436865550e+00;
  value[6].exponentP[0] = 0.3049723950e+02;
  value[6].exponentP[1] = 0.6036199601e+01;
  value[6].exponentP[2] = 0.1876046337e+01;
  value[6].exponentP[3] = 0.7217826470e+00;
  value[6].exponentP[4] = 0.3134706954e+00;
  value[6].exponentP[5] = 0.1436865550e+00;
  value[6].coeffS[0] = -0.1325278809e-01;
  value[6].coeffS[1] = -0.4699171014e-01;
  value[6].coeffS[2] = -0.3378537151e-01;
  value[6].coeffS[3] =  0.2502417861e+00;
  value[6].coeffS[4] =  0.5951172526e+00;
  value[6].coeffS[5] =  0.2407061763e+00;
  value[6].coeffP[0] =  0.3759696623e-02;
  value[6].coeffP[1] =  0.3767936984e-01;
  value[6].coeffP[2] =  0.1738967435e+00;
  value[6].coeffP[3] =  0.4180364347e+00;
  value[6].coeffP[4] =  0.4258595477e+00;
  value[6].coeffP[5] =  0.1017082955e+00;
  /***************************************************************************************/ 
  //
  // Stolen from Sparrow
  //
  /***************************************************************************************/ 
   value[6].exponentS[0] =   8.846144274163005e+01;  
   value[6].exponentS[1] =   1.622292919221501e+01;  
   value[6].exponentS[2] =   4.558993095148346e+00;  
   value[6].exponentS[3] =   6.519461344814428e-01;  
   value[6].exponentS[4] =   2.958932522507792e-01;  
   value[6].exponentS[5] =   1.411098199523801e-01;  
  
      value[6].coeffS[0] =  -8.534093658582062e-02;
      value[6].coeffS[1] =  -1.190837759635380e-01;
      value[6].coeffS[2] =  -1.145234971659361e-01;
      value[6].coeffS[3] =   1.730334274890658e-01;
      value[6].coeffS[4] =   1.607236215755646e-01;
      value[6].coeffS[5] =   2.810826497961667e-02;
  
    value[6].exponentP[0] =  1.875086664726173e+01;  
    value[6].exponentP[1] =  4.889844711496806e+00;  
    value[6].exponentP[2] =  1.749633028685196e+00;  
    value[6].exponentP[3] =  7.313800499385337e-01;  
    value[6].exponentP[4] =  3.344367808801327e-01;  
    value[6].exponentP[5] =  1.581099099775126e-01;  
                          
       value[6].coeffP[0] =  4.407314493338078e-01;
       value[6].coeffP[1] =  5.331737162464377e-01;
       value[6].coeffP[2] =  5.445171182695406e-01;
       value[6].coeffP[3] =  3.904445162839797e-01;
       value[6].coeffP[4] =  1.454563993779586e-01;
       value[6].coeffP[5] =  1.494840672526637e-02;
  /***************************************************************************************/ 
  //
  // End Stolen from Sparrow
  //
  /***************************************************************************************/ 
/***************************************************************************************/ 
// Nytrogen 2s,2p
  value[7].exponentS[0] =   0.3919880787e+02; 
  value[7].exponentS[1] =   0.7758467071e+01; 
  value[7].exponentS[2] =   0.2411325783e+01; 
  value[7].exponentS[3] =   0.9277239437e+00; 
  value[7].exponentS[4] =   0.4029111410e+00; 
  value[7].exponentS[5] =   0.1846836552e+00; 
  value[7].coeffS[0]    =  -0.1325278809e-01;  
  value[7].coeffS[1]    =  -0.4699171014e-01;  
  value[7].coeffS[2]    =  -0.3378537151e-01;  
  value[7].coeffS[3]    =   0.2502417861e+00;  
  value[7].coeffS[4]    =   0.5951172526e+00;  
  value[7].coeffS[5]    =   0.2407061763e+00;  
  value[7].exponentP[0] =   0.3919880787e+02; 
  value[7].exponentP[1] =   0.7758467071e+01; 
  value[7].exponentP[2] =   0.2411325783e+01; 
  value[7].exponentP[3] =   0.9277239437e+00; 
  value[7].exponentP[4] =   0.4029111410e+00; 
  value[7].exponentP[5] =   0.1846836552e+00; 
  value[7].coeffP[0]    =   0.3759696623e-02;
  value[7].coeffP[1]    =   0.3767936984e-01;
  value[7].coeffP[2]    =   0.1738967435e+00;
  value[7].coeffP[3]    =   0.4180364347e+00;
  value[7].coeffP[4]    =   0.4258595477e+00;
  value[7].coeffP[5]    =   0.1017082955e+00;
  /***************************************************************************************/ 
  //
  // Stolen from Sparrow
  //
  /***************************************************************************************/ 
   value[7].exponentS[0] =   1.408553999528384e+02;  
   value[7].exponentS[1] =   2.583144824406823e+01;  
   value[7].exponentS[2] =   7.259194242116380e+00;  
   value[7].exponentS[3] =   1.038080893483724e+00;  
   value[7].exponentS[4] =   4.711449542017903e-01;  
   value[7].exponentS[5] =   2.246863662931398e-01;  
  
      value[7].coeffS[0] =  -1.209684224448265e-01;
      value[7].coeffS[1] =  -1.687979660569579e-01;
      value[7].coeffS[2] =  -1.623338967119989e-01;
      value[7].coeffS[3] =   2.452701082384343e-01;
      value[7].coeffS[4] =   2.278212980714544e-01;
      value[7].coeffS[5] =   3.984268990094886e-02;
  
    value[7].exponentP[0] =  2.985663289232627e+01;  
    value[7].exponentP[1] =  7.786002705798312e+00;  
    value[7].exponentP[2] =  2.785905953918744e+00;  
    value[7].exponentP[3] =  1.164561940873007e+00;  
    value[7].exponentP[4] =  5.325170500259336e-01;  
    value[7].exponentP[5] =  2.517552723103985e-01;  
                          
       value[7].coeffP[0] =  7.883124678886093e-01;
       value[7].coeffP[1] =  9.536589428843584e-01;
       value[7].coeffP[2] =  9.739482715823707e-01;
       value[7].coeffP[3] =  6.983669549124415e-01;
       value[7].coeffP[4] =  2.601699818271488e-01;
       value[7].coeffP[5] =  2.673740531657013e-02;
  /***************************************************************************************/ 
  //
  // End Stolen from Sparrow
  //
  /***************************************************************************************/ 
/***************************************************************************************/ 
// Oxygen 2s,2p
  value[8].exponentS[0] = 0.5218776196e+02;   
  value[8].exponentS[1] = 0.1032932006e+02;   
  value[8].exponentS[2] = 0.3210344977e+01;   
  value[8].exponentS[3] = 0.1235135428e+01;   
  value[8].exponentS[4] = 0.5364201581e+00;   
  value[8].exponentS[5] = 0.2458806060e+00;   
  value[8].coeffS[0] =  -0.1325278809e-01;  
  value[8].coeffS[1] =  -0.4699171014e-01;  
  value[8].coeffS[2] =  -0.3378537151e-01;  
  value[8].coeffS[3] =   0.2502417861e+00;  
  value[8].coeffS[4] =   0.5951172526e+00;  
  value[8].coeffS[5] =   0.2407061763e+00;  
  value[8].exponentP[0] = 0.5218776196e+02;   
  value[8].exponentP[1] = 0.1032932006e+02;   
  value[8].exponentP[2] = 0.3210344977e+01;   
  value[8].exponentP[3] = 0.1235135428e+01;   
  value[8].exponentP[4] = 0.5364201581e+00;   
  value[8].exponentP[5] = 0.2458806060e+00;   
  value[8].coeffP[0] =   0.3759696623e-02;
  value[8].coeffP[1] =   0.3767936984e-01;
  value[8].coeffP[2] =   0.1738967435e+00;
  value[8].coeffP[3] =   0.4180364347e+00;
  value[8].coeffP[4] =   0.4258595477e+00;
  value[8].coeffP[5] =   0.1017082955e+00;
  
  /***************************************************************************************/ 
  //
  // Stolen from Sparrow
  //
  /***************************************************************************************/ 
   value[8].exponentS[0] =    2.018091738330405e+02; 
   value[8].exponentS[1] =    3.700975064350955e+01; 
   value[8].exponentS[2] =    1.040053837613310e+01; 
   value[8].exponentS[3] =    1.487300079059505e+00; 
   value[8].exponentS[4] =    6.750282487920551e-01; 
   value[8].exponentS[5] =    3.219171573709536e-01; 
                              
      value[8].coeffS[0] =   -1.584154624760238e-01;
      value[8].coeffS[1] =   -2.210511414259479e-01;
      value[8].coeffS[2] =   -2.125860518260085e-01;
      value[8].coeffS[3] =    3.211960348235300e-01;
      value[8].coeffS[4] =    2.983457630220672e-01;
      value[8].coeffS[5] =    5.217641116074215e-02;
  
   value[8].exponentP[0] =     4.277679392805792e+01;  
   value[8].exponentP[1] =     1.115531796470053e+01;  
   value[8].exponentP[2] =     3.991479056714438e+00;  
   value[8].exponentP[3] =     1.668514542173559e+00;  
   value[8].exponentP[4] =     7.629585088943974e-01;  
   value[8].exponentP[5] =     3.606998633356254e-01; 
                          
      value[8].coeffP[0] =     1.235684643746845e+00;
      value[8].coeffP[1] =     1.494866260646997e+00;
      value[8].coeffP[2] =     1.526669908217376e+00;
      value[8].coeffP[3] =     1.094694498739664e+00;
      value[8].coeffP[4] =     4.078180472887439e-01;
      value[8].coeffP[5] =     4.191104734371579e-02;
  
  /***************************************************************************************/ 
  //
  // End Stolen from Sparrow
  //
  /***************************************************************************************/ 
}
/***************************************************************************************/
void STO_6G::To_device(){
 // printf ("STO_6G To_device this = %p\n",this);
  #pragma acc enter data copyin(this[0:1])
  #pragma acc enter data copyin(this->value[0:countAtomValues_])
  for (unsigned int i=0;i<countAtomValues_;++i) {
    #pragma acc enter data copyin(this->value[i].exponentS[0:6],\
        this->value[i].exponentP[0:6],\
        this->value[i].coeffS[0:6],\
        this->value[i].coeffP[0:6])
  }
  //printf ("STO_6G To_device this = %p\n",this);
}
/***************************************************************************************/ 
void STO_6G::From_device(){
  for (unsigned int i=0;i<countAtomValues_;++i) {
    #pragma acc exit data copyout(this->value[i].exponentS[0:6],\
        this->value[i].exponentP[0:6],\
        this->value[i].coeffS[0:6],\
        this->value[i].coeffP[0:6])
  }
  #pragma acc exit data copyout(this[0:1])
}
#endif // _STO-6G_CPP_
