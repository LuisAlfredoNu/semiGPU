#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "mathtools.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for MathTool Class " << endl;
	cout << "********************************************************" << endl << endl;
	
   MathTool mathtool;

   // Carteian coordinates
	double x1=1.5, y1=2.0, z1=0.5;
   // Spherical coordinates in Rad
	double r=2.2, theta=2.4, phi=2.6;

   vector<double> xyz_coor = {x1,y1,z1};
   vector<double> spherical_coor = {r,theta,phi};

   vector<double> xyz2spherical = mathtool.changeCartesian2Spherical(xyz_coor);
   vector<double> spherical2xyz = mathtool.changeSpherical2Cartesian(spherical_coor);

   cout << "Change coordinates system " << endl << endl;
   cout << "Cartesian to Spherical " << endl;
   for(int i=0;i<3;i++){
      cout << xyz_coor[i] << " -> " << xyz2spherical[i] << endl;
   }
   cout << "Spherical to Cartesian " << endl;
   for(int i=0;i<3;i++){
      cout << spherical_coor[i] << " -> " << spherical2xyz[i] << endl;
   }
	return EXIT_SUCCESS;
}


