/*
 * List of STO-6G basis
*/
#ifndef _STO_6G_H_
#define _STO_6G_H_

/***************************************************************************************/
class STO_6G{
/***************************************************************************************/ 
 public:
  STO_6G();

  struct STOdata{
    double exponent[6];
    double coeffS[6];
    double coeffP[6];
  };

  struct STOdata value[36];
};


#endif // _STO-6G_H_
