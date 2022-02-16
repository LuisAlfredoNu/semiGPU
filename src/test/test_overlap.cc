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

#include "basematrix.h"
#include "overlap.h"
int main (int argc, char *argv[]){
  cout << endl;
  cout << "********************************************************" << endl;
  cout << " Testing Class: Overlap " << endl;
  cout << "********************************************************" << endl << endl;
/***************************************************************************************/ 
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

  /**/ 
  vector<Atom> molecule (6,Atom());
  for (int i=0;i<6;++i) {
    molecule[i].setAtomNumber(1);
  }
  double coorA[3] = { 0.00000, 0.00000, 0.00000};  
  double coorB[3] = { 0.98320, 0.00000, 0.00000};
  double coorC[3] = { 1.47480,-0.85148, 0.00000};
  double coorD[3] = { 0.98320,-1.70290, 0.00000};
  double coorF[3] = { 0.00000,-1.70290, 0.00000};
  double coorG[3] = {-0.49160,-0.85148, 0.00000};

  molecule[0].setCoordinates(coorA);
  molecule[1].setCoordinates(coorB);
  molecule[2].setCoordinates(coorC);
  molecule[3].setCoordinates(coorD);
  molecule[4].setCoordinates(coorF);
  molecule[5].setCoordinates(coorG);

  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);

  Overlap overlapA(infoAOs);

  if (overlapA.Alloc4Matrix("overlapMatrix")){
    cout << "Correct Alloc: overlapMatrix" << endl;
  }else{
    cout << "Bad Alloc: overlapMatrix " << endl;
  }

  overlapA.ComputeMatrix();

  cout << "Overlap Matrix done" << endl;
  cout << "HexaHidrogen" << endl;
  ScreenUtils::PrintMatrixNxNSymmetric(infoAOs.orbital.size(),overlapA.matrixHold_);

  overlapValue = overlap.ComputeOverlap(infoAOs.orbital[0],infoAOs.orbital[0]);

  cout << "overlapValue = " << overlapValue << endl;
  size_t As = 0;
  size_t Bs = 0;
  overlapValue = overlap.ComputeOverlap(infoAOs.orbital[As],infoAOs.orbital[Bs]);

  cout << "overlapValue = " << overlapValue << endl;

  vector<Atom> moleculeA (2,Atom());
  double coorAAA[3] = {-1.00000,-0.00000,-0.00000};  
  double coorBBB[3] = { 0.00000, 0.00000, 0.50000};

  moleculeA[0].setAtomNumber(6);
  moleculeA[0].setCoordinates(coorAAA[0],coorAAA[1],coorAAA[2]);
  moleculeA[1].setAtomNumber(6);
  moleculeA[1].setCoordinates(coorBBB[0],coorBBB[1],coorBBB[2]);


  ListAtomicOrbitals infoAOsA;
  infoAOsA.SetOrbitals(moleculeA);

  Overlap overlapB(infoAOsA);

  if (overlapB.Alloc4Matrix("overlapMatrix")){
    cout << "Correct Alloc: overlapMatrix" << endl;
  }else{
    cout << "Bad Alloc: overlapMatrix " << endl;
  }

  overlapB.ComputeMatrix();

  cout << "Overlap Matrix done" << endl;
  cout << "1C -> {"<<coorAAA[0]<<","<<coorAAA[1]<<","<<coorAAA[2]<<"}" << endl;
  cout << "2C -> {"<<coorBBB[0]<<","<<coorBBB[1]<<","<<coorBBB[2]<<"}" << endl;
  ScreenUtils::PrintMatrixNxNSymmetric(infoAOsA.orbital.size(),overlapB.matrixHold_);
  /**/
  return EXIT_SUCCESS;
}


