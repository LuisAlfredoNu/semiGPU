#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <iomanip>
using std::setw;
using std::setprecision;
#include <string>
using std::string;
#include <vector>
using std::vector;

#include "atom.h"
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
  int angularMomA[3] = {1,0,0};
  orbitalA.SetElement( 6);
  orbitalA.SetCoordinates(coorA);
  orbitalA.SetAngularMomentum(angularMomA);
  orbitalA.SetIndexAtom(0);

  int angularMomB[3] = {1,0,0};
  orbitalB.SetElement( 6);
  orbitalB.SetCoordinates(coorA);
  orbitalB.SetAngularMomentum(angularMomB);
  orbitalB.SetIndexAtom(0);

  double coorC[3] = {1.0,1.0,1.5};
  int angularMomC[3] = {0,0,0};
  orbitalC.SetElement( 6);
  orbitalC.SetCoordinates(coorC);
  orbitalC.SetAngularMomentum(angularMomC);
  orbitalC.SetIndexAtom(1);

  int angularMomD[3] = {1,0,0};
  orbitalD.SetElement( 6);
  orbitalD.SetCoordinates(coorC);
  orbitalD.SetAngularMomentum(angularMomD);
  orbitalD.SetIndexAtom(1);

  cout << "orbitalA.angularMomentumInt = " << orbitalA.angularMomentumInt << endl;
  cout << "orbitalB.angularMomentumInt = " << orbitalB.angularMomentumInt << endl;
  cout << "orbitalC.angularMomentumInt = " << orbitalC.angularMomentumInt << endl;
  cout << "orbitalD.angularMomentumInt = " << orbitalD.angularMomentumInt << endl;
  cout << endl;
  // Init MNDO parameters
  MNDOparameter MNDOpara;

  // Init TwoCenterIntegral
  TwoCenterIntegral twoCIntegral(MNDOpara);

  double result =  twoCIntegral.ComputeTwoCenterIntegral(orbitalA,orbitalB,orbitalC,orbitalD);
  cout << "result = " << result << endl;

  vector<Atom> molecule (2,Atom());

  molecule[0].setCoordinates(coorA[0],coorA[1],coorA[2]);
  molecule[0].setAtomNumber(6);
  molecule[1].setCoordinates(coorC[0],coorC[1],coorC[2]);
  molecule[1].setAtomNumber(6);

  ListAtomicOrbitals AOs;
  AOs.SetOrbitals(molecule);

  double**** all2CenterIntegral = NULL;
  if (TwoCenterIntegral::Alloc4AllTwoCenterIntegral(molecule,all2CenterIntegral)){
    cout << "Correct Alloc" << endl;
  }else{
    cout << "Bad Alloc" << endl;
  }

  twoCIntegral.ComputeAllTwoCenterIntegral(AOs.orbital,all2CenterIntegral);

  cout << std::fixed << setprecision(4);
  for (int i=0;i<2;i++) {
    for (int j=0;j<=i;j++) {
      for (int k=0;k<10;k++) {
        for (int l=0;l<10;l++) {
          cout << setw(12) << all2CenterIntegral[i][j][k][l];
        }
        cout << endl;
      }
        cout << endl;
    }
  }

	return EXIT_SUCCESS;
}


