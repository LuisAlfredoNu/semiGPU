#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;
#include <vector>
using std::vector;

#include "atomicOrbitals.h"
#include "MNDO_parameters.h"

#include "twocenterintegral.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for Class TwoCenterIntegral " << endl;
	cout << "********************************************************" << endl << endl;


  // Make singles orbitals

  AtomicOrbital orbitalA,orbitalB,orbitalC,orbitalD;

  double coorA[3] = {0.0,0.0,0.0};
  int angularMomA[3] = {0,0,0};
  orbitalA.SetElement( 6);
  orbitalA.SetCoordinates(coorA);
  orbitalA.SetAngularMomentum(angularMomA);
  orbitalA.SetNAtom(0);

  int angularMomB[3] = {0,0,0};
  orbitalB.SetElement( 6);
  orbitalB.SetCoordinates(coorA);
  orbitalB.SetAngularMomentum(angularMomB);
  orbitalB.SetNAtom(0);

  double coorC[3] = {1.0,1.0,1.5};
  int angularMomC[3] = {0,0,0};
  orbitalC.SetElement( 6);
  orbitalC.SetCoordinates(coorC);
  orbitalC.SetAngularMomentum(angularMomC);
  orbitalC.SetNAtom(1);

  int angularMomD[3] = {0,0,0};
  orbitalD.SetElement( 6);
  orbitalD.SetCoordinates(coorC);
  orbitalD.SetAngularMomentum(angularMomD);
  orbitalD.SetNAtom(1);

  // Init MNDO parameters
  MNDOparameter MNDOpara;

  // Init TwoCenterIntegral
  TwoCenterIntegral twoCIntegral(MNDOpara);

  double result =  twoCIntegral.ComputeTwoCenterIntegral(orbitalA,orbitalB,orbitalC,orbitalD);
  cout << "result = " << result << endl;


	return EXIT_SUCCESS;
}


