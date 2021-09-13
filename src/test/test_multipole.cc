#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;
#include <vector>
using std::vector;

#include "MNDO_parameters.h"
#include "mymath.h"

#include "multipole.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for Class Multipole " << endl;
	cout << "********************************************************" << endl << endl;


  // Make singles orbitals

  AtomicOrbital orbitalA,orbitalB;

  double coorAA[3] = {0.0,0.0,0.0};
  orbitalA.SetElement( 6);
  orbitalA.SetCoordinates(coorAA);

  double coorBB[3] = {0.0,0.0,1.0};
  orbitalB.SetElement( 8);
  orbitalB.SetCoordinates(coorBB);

  // Init MNDO parameters
  MNDOparameter MNDOpara;
  Multipole multipole(MNDOpara);

  double atomDistance = distancePointsV3(orbitalA.coordinates, orbitalB.coordinates);
  cout << "atomDistance = " << atomDistance << endl;

  double interaction;

  interaction = multipole.Interaction_qq(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction qq = " << interaction << endl;

  interaction = multipole.Interaction_qUz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction qUz = " << interaction << endl;

  interaction = multipole.Interaction_Uzq(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction Uzq = " << interaction << endl;

  interaction = multipole.Interaction_qQpipi(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction qQpipi = " << interaction << endl;

  interaction = multipole.Interaction_Qpipiq(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction Qpipiq = " << interaction << endl;

  interaction = multipole.Interaction_qQzz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction qQzz = " << interaction << endl;

  interaction = multipole.Interaction_Qzzq(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction Qzzq = " << interaction << endl;

  interaction = multipole.Interaction_UpiUpi(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction UpiUpi = " << interaction << endl;

  interaction = multipole.Interaction_UzUz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction UzUz = " << interaction << endl;

  cout << "******************************" << endl;
  interaction = multipole.Interaction_UpiQpiz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction UpiQpiz = " << interaction << endl;

  interaction = multipole.Interaction_QpizUpi(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QpizUpi = " << interaction << endl;

  interaction = multipole.Interaction_UzQpipi(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction UzQpipi = " << interaction << endl;

  interaction = multipole.Interaction_QpipiUz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QpipiUz = " << interaction << endl;

  interaction = multipole.Interaction_UzQz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction UzQz = " << interaction << endl;

  interaction = multipole.Interaction_QzUz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QzUz = " << interaction << endl;

  interaction = multipole.Interaction_QpipiQpipi(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QpipiQpipi = " << interaction << endl;

  interaction = multipole.Interaction_QxxQyy(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QxQy = " << interaction << endl;

  interaction = multipole.Interaction_QpipiQzz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QpipiQzz = " << interaction << endl;

  interaction = multipole.Interaction_QzzQpipi(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QzzQpipi = " << interaction << endl;

  cout << "******************************" << endl;
  
  interaction = multipole.Interaction_QzzQzz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QzzQzz = " << interaction << endl;
  
  interaction = multipole.Interaction_QpizQpiz(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QpizQpiz = " << interaction << endl;

  interaction = multipole.Interaction_QxyQxy(atomDistance,orbitalA,orbitalB);
  cout << "multipole interaction QxyQxy = " << interaction << endl;

	return EXIT_SUCCESS;
}


