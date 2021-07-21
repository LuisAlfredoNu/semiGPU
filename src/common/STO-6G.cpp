#ifndef _STO_6G_CPP_
#define _STO_6G_CPP_

#include "STO-6G.h"
#include <iostream>
using std::cout;

STO_6G::STO_6G(){

// Hidrogen 1s
value[1].exponent[0] = 0.3552322122e+02;       
value[1].exponent[1] = 0.6513143725e+01;       
value[1].exponent[2] = 0.1822142904e+01;       
value[1].exponent[3] = 0.6259552659e+00;       
value[1].exponent[4] = 0.2430767471e+00;       
value[1].exponent[5] = 0.1001124280e+00;      

value[1].coeffS[0] = 0.9163596281e-02;
value[1].coeffS[1] = 0.4936149294e-01;
value[1].coeffS[2] = 0.1685383049e+00;
value[1].coeffS[3] = 0.3705627997e+00;
value[1].coeffS[4] = 0.4164915298e+00;
value[1].coeffS[5] = 0.1303340841e+00;

/***************************************************************************************/ 
// Carbon 2s, 2p
value[6].exponent[0] = 0.3049723950e+02;
value[6].exponent[1] = 0.6036199601e+01;
value[6].exponent[2] = 0.1876046337e+01;
value[6].exponent[3] = 0.7217826470e+00;
value[6].exponent[4] = 0.3134706954e+00;
value[6].exponent[5] = 0.1436865550e+00;
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

};
#endif // _STO-6G_CPP_
