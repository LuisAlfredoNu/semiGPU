/*grep -r "STO-6G"
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

#include "atomicOrbitals.h"
/***************************************************************************************/
class Overlap{
 public:
  Overlap();
/***************************************************************************************/ 
  // Variables
  // Start STO-6G data
  STO_6G basisSTO;
/***************************************************************************************/ 
  // Methods 
  bool GetOverlapMatrix(const vector<AtomicOrbital> infoOrbitals,double* &overlapMatrix);

  double ComputeOverlap(const AtomicOrbital,const AtomicOrbital);
 private:

};
#endif // _OVERLAP_H_


