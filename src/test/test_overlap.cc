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
#include "screenutils.h"

#include "overlap.h"
int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing Class: Overlap " << endl;
	cout << "********************************************************" << endl << endl;

  // Make singles orbitals

  AtomicOrbital orbitalA,orbitalB;

  double coorAA[3] = {0.0,0.0,0.0};
  orbitalA.SetElement( 1);
  orbitalA.SetCoordinates(coorAA);

  double coorBB[3] = {0.0,0.0,1.0};
  orbitalB.SetElement( 1);
  orbitalB.SetCoordinates(coorBB);

  Overlap overlap;
  double overlapValue = overlap.ComputeOverlap(orbitalA,orbitalB);

  cout << " Overlap between H-1s H-1s with 1 Angstrom distance" << endl;
  cout << " Expected value = 0.496735 " << endl; 
  cout << "  Compute value = " << overlapValue << endl << endl;

  double coorAAAA[3] = {0.0,0.0,0.0};
  orbitalA.SetCoordinates(coorAAAA);

  double coorBBBB[3] = {1.0,1.0,1.0};
  orbitalB.SetCoordinates(coorBBBB);
  
  orbitalA.SetElement( 6);
  int angularMomentumA[3] = {0,0,0};
  orbitalA.SetAngularMomentum(angularMomentumA);
  
  orbitalB.SetElement( 6);
  int angularMomentumB[3] = {1,0,0};
  orbitalB.SetAngularMomentum(angularMomentumB);
  
  overlapValue = overlap.ComputeOverlap(orbitalA,orbitalB);

  cout << " Overlap between C-2s C-2p_y with 1 Angstrom distance" << endl;
  cout << " Expected value = -0.152276 " << endl; 
  cout << "  Compute value = " << overlapValue << endl << endl;
  
  int angularMomentumAA[3] = {1,0,0};
  orbitalA.SetAngularMomentum(angularMomentumAA);
  
  int angularMomentumBB[3] = {1,0,0};
  orbitalB.SetAngularMomentum(angularMomentumBB);
  
  overlapValue = overlap.ComputeOverlap(orbitalA,orbitalB);

  cout << " Overlap between C-2p_x C-2p_x with 1 Angstrom distance" << endl;
  cout << " Expected value = -0.0184144 " << endl; 
  cout << "  Compute value = " << overlapValue << endl << endl;
  
  int angularMomentumAAA[3] = {1,0,0};
  orbitalA.SetAngularMomentum(angularMomentumAAA);
  
  int angularMomentumBBB[3] = {0,0,1};
  orbitalB.SetAngularMomentum(angularMomentumBBB);
  
  overlapValue = overlap.ComputeOverlap(orbitalA,orbitalB);

  cout << " Overlap between C-2p_x C-2p_z with 1 Angstrom distance" << endl;
  cout << " Expected value = -0.130419 " << endl; 
  cout << "  Compute value = " << overlapValue << endl << endl;
  
  vector<AtomicOrbital> infoOrbitals (6,AtomicOrbital());

  double coorA[3] = { 0.00000, 0.00000, 0.00000};  
  infoOrbitals[0].SetElement(1);
  infoOrbitals[0].SetCoordinates(coorA);
  double coorB[3] = { 0.98320, 0.00000, 0.00000};
  infoOrbitals[1].SetElement(1);
  infoOrbitals[1].SetCoordinates(coorB);
  double coorC[3] = { 1.47480,-0.85148, 0.00000};
  infoOrbitals[2].SetElement(1);
  infoOrbitals[2].SetCoordinates(coorC);
  double coorD[3] = { 0.98320,-1.70290, 0.00000};
  infoOrbitals[3].SetElement(1);
  infoOrbitals[3].SetCoordinates(coorD);
  double coorF[3] = { 0.00000,-1.70290, 0.00000};
  infoOrbitals[4].SetElement(1);
  infoOrbitals[4].SetCoordinates(coorF);
  double coorG[3] = {-0.49160,-0.85148, 0.00000};
  infoOrbitals[5].SetElement(1);
  infoOrbitals[5].SetCoordinates(coorG);
  
  double* overlapMatrix;
  bool doneOverlap = overlap.GetOverlapMatrix(infoOrbitals,overlapMatrix);

  if (doneOverlap) {
    cout << "Overlap Matrix done" << endl;
    cout << "HexaHidrogen" << endl;
    ScreenUtils::PrintMatrixNxNSymmetric(infoOrbitals.size(),overlapMatrix);
  } else {
    string err = "Fail overlap test";
    ScreenUtils::DisplayErrorMessage(err);
  } 

  vector<AtomicOrbital> infoOrbitalsB (8,AtomicOrbital());

  double coorAAA[3] = {-0.30000,-0.30000,-0.30000};  
  double coorBBB[3] = { 0.30000, 0.20000, 0.50000};
  infoOrbitalsB[0].SetElement(6);
  infoOrbitalsB[0].SetCoordinates(coorAAA);
  angularMomentumA[0] = 0;
  angularMomentumA[1] = 0;
  angularMomentumA[2] = 0;
  infoOrbitalsB[0].SetAngularMomentum(angularMomentumA);

  infoOrbitalsB[1].SetElement(6);
  infoOrbitalsB[1].SetCoordinates(coorAAA);
  angularMomentumA[0] = 1;
  angularMomentumA[1] = 0;
  angularMomentumA[2] = 0;
  infoOrbitalsB[1].SetAngularMomentum(angularMomentumA);

  infoOrbitalsB[2].SetElement(6);
  infoOrbitalsB[2].SetCoordinates(coorAAA);
  angularMomentumA[0] = 0;
  angularMomentumA[1] = 1;
  angularMomentumA[2] = 0;
  infoOrbitalsB[2].SetAngularMomentum(angularMomentumA);

  infoOrbitalsB[3].SetElement(6);
  infoOrbitalsB[3].SetCoordinates(coorAAA);
  angularMomentumA[0] = 0;
  angularMomentumA[1] = 0;
  angularMomentumA[2] = 1;
  infoOrbitalsB[3].SetAngularMomentum(angularMomentumA);

  infoOrbitalsB[4].SetElement(6);
  infoOrbitalsB[4].SetCoordinates(coorBBB);
  angularMomentumA[0] = 0;
  angularMomentumA[1] = 0;
  angularMomentumA[2] = 0;
  infoOrbitalsB[4].SetAngularMomentum(angularMomentumA);

  infoOrbitalsB[5].SetElement(6);
  infoOrbitalsB[5].SetCoordinates(coorBBB);
  angularMomentumA[0] = 1;
  angularMomentumA[1] = 0;
  angularMomentumA[2] = 0;
  infoOrbitalsB[5].SetAngularMomentum(angularMomentumA);

  infoOrbitalsB[6].SetElement(6);
  infoOrbitalsB[6].SetCoordinates(coorBBB);
  angularMomentumA[0] = 0;
  angularMomentumA[1] = 1;
  angularMomentumA[2] = 0;
  infoOrbitalsB[6].SetAngularMomentum(angularMomentumA);

  infoOrbitalsB[7].SetElement(6);
  infoOrbitalsB[7].SetCoordinates(coorBBB);
  angularMomentumA[0] = 0;
  angularMomentumA[1] = 0;
  angularMomentumA[2] = 1;
  infoOrbitalsB[7].SetAngularMomentum(angularMomentumA);
  
  double* overlapMatrixB;
   doneOverlap = overlap.GetOverlapMatrix(infoOrbitalsB,overlapMatrixB);

  if (doneOverlap) {
    cout << "Overlap Matrix done" << endl;
    cout << "1C -> {0.0,0.0,0.0}" << endl;
    cout << "2C -> {1.0,1.0,1.0}" << endl;
    ScreenUtils::PrintMatrixNxNSymmetric(infoOrbitalsB.size(),overlapMatrixB);
  } else {
    string err = "Fail overlap test";
    ScreenUtils::DisplayErrorMessage(err);
  }
	return EXIT_SUCCESS;
}


