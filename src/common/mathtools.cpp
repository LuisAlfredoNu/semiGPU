/*
 * Class with mathematics tools 
 */

#ifndef _MATHTOOLS_CPP_
#define _MATHTOOLS_CPP_

#include <vector>
using std::vector;
#include <math.h> 

#include "mathtools.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
MathTool::MathTool(){}
/***************************************************************************************/ 
vector<double> MathTool::changeCartesian2Spherical(vector<double> xyz_coor){
   vector<double> spherical_coor (3,0.0);

   // r coordinate
   spherical_coor[0] = sqrt(xyz_coor[0]*xyz_coor[0] + xyz_coor[1]*xyz_coor[1] +xyz_coor[2]*xyz_coor[2]);
   // theta coordinate
   spherical_coor[1] = acos(xyz_coor[2]/spherical_coor[0]); 
   // phi coordinate
   spherical_coor[2] = atan2(xyz_coor[1],xyz_coor[0]);


   return spherical_coor;
}
/***************************************************************************************/ 

vector<double> MathTool::changeSpherical2Cartesian(vector<double> spherical_coor){
   vector<double> xyz_coor (3,0.0);

   // x coordinate
   xyz_coor[0] = spherical_coor[0] * sin(spherical_coor[1]) * cos(spherical_coor[2]) ;
   // y coordinate
   xyz_coor[1] = spherical_coor[0] * sin(spherical_coor[1]) * sin(spherical_coor[2]) ;
   // z coordinate
   xyz_coor[2] = spherical_coor[0] * cos(spherical_coor[1]) ;

   return xyz_coor;
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _MATHTOOLS_CPP_ 
