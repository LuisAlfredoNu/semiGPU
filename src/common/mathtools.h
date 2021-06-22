/*
 * Class with mathematics tools 
 */

#ifndef _MATHTOOLS_H_
#define _MATHTOOLS_H_

#include <vector>
using std::vector;
/***************************************************************************************/ 

class MathTool{
/***************************************************************************************/ 
public:
   MathTool();
/***************************************************************************************/ 
   vector<double> changeCartesian2Spherical(vector<double>);
   vector<double> changeSpherical2Cartesian(vector<double>);
/***************************************************************************************/ 
private:
};

#endif // _MATHTOOLS_H_ 


