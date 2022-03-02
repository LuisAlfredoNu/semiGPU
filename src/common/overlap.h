/*
 * Class for compute the overlap between 2 orbitals
 * -- Develop the more basic methed: Boys
 * -- Numerical method Ref: https://doi.org/10.1016/0097-8485(78)80005-9 Ec. 29 - 33 
 * -- Implementation Ref: https://www.mathematica-journal.com/2012/02/16/evaluation-of-gaussian-molecular-integrals/
*/

#ifndef _OVERLAP_H_
#define _OVERLAP_H_

#include <vector>
using std::vector;

#include "STO-6G.h"

#include "basematrix.h"
#include "atomicOrbitals.h"
/***************************************************************************************/
class Overlap: public BaseMatrix{
 public:
  Overlap();
  // Constructor only to create a overlap matrix
  Overlap(const ListAtomicOrbitals &);
/***************************************************************************************/ 
  // Variables
  // Start STO-6G data
  STO_6G* basisSTO;

/***************************************************************************************/ 
  // Methods 
  double ComputeElementMatrix(const size_t &i,const size_t &j);
  double ComputeOverlap(const AtomicOrbital&,const AtomicOrbital&);
  double ComputeOverlap_Boys(const AtomicOrbital&,const AtomicOrbital&);

 private:
/***************************************************************************************/ 
  // Variables
  const ListAtomicOrbitals *infoAOs_;
  int tmp_count;
  double *xk_, *wk_;
/***************************************************************************************/ 
  // Methods 
  // Overlap McMurchieDavidson Method
  double McMurchieCoef(const int &,const int &,const int &,const double &,const double &,\
         const double &);
  double OverlapMcMurchie(const int &,const int &,const double &,const double &,const double &);
  double NormalizationConst(const double &,const int (&)[3]);
  // Overlap Numerical Integral
  double OverlapNumericalIntegral(const int &,const int &,const double &,const double &,const double &);
  // Overlap Boys Method
  void Overlap_SS(double&,const double&,const double&,const double&);
  void Overlap_PS(double&,const double&,const double&,const double&,const double&);
  void Overlap_SP(double&,const double&,const double&,const double&,const double&);
  void Overlap_PxPy(double&,const double&,const double&,const double&,const double&);
  void Overlap_PxPx(double&,const double&,const double&,const double&);

};
#endif // _OVERLAP_H_


