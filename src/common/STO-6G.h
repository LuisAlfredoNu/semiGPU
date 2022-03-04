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
  ~STO_6G();

  struct STOdata{
    double exponentS[6];
    double exponentP[6];
    double coeffS[6];
    double coeffP[6];
  };

  unsigned int countAtomValues_;
  struct STOdata* value;
  private:
/***************************************************************************************/
// Methods
/***************************************************************************************/
   void SetArraysSizes();
   void SetParameterValues();
   void To_device();
   void From_device();
};


#endif // _STO-6G_H_
