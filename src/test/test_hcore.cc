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
#include "screenutils.h"

#include "hcore.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for Class Hcore " << endl;
	cout << "********************************************************" << endl << endl;

  // Init MNDO parameters
  MNDOparameter MNDOpara;

  // Init molecule
  vector<Atom> molecule (2,Atom());

  double coorA[3] = {-1.0,-2.0,-3.0};
  double coorC[3] = {1.0,1.0,1.5};
  molecule[0].setCoordinates(coorA[0],coorA[1],coorA[2]);
  molecule[0].setAtomNumber(6);
  molecule[1].setCoordinates(coorC[0],coorC[1],coorC[2]);
  molecule[1].setAtomNumber(6);

  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);

  // Alloc for all 2 center integrals 
  double**** all2CenterIntegral = NULL;
  if (TwoCenterIntegral::Alloc4AllTwoCenterIntegral(molecule,all2CenterIntegral)){
    cout << "Correct Alloc: all2CenterIntegral" << endl;
  }else{
    cout << "Bad Alloc: all2CenterIntegral " << endl;
  }
  // Init TwoCenterIntegral
  TwoCenterIntegral twoCIntegral(MNDOpara);
  twoCIntegral.ComputeAllTwoCenterIntegral(infoAOs,all2CenterIntegral);

  // Alloc Hcore Matrix
  double* hcoreMatrix = NULL;
  if (Hcore::Alloc4HcoreMatrix(infoAOs.orbital.size(),hcoreMatrix)){
    cout << "Correct Alloc: hcoreMatrix" << endl;
  }else{
    cout << "Bad Alloc: hcoreMatrix " << endl;
  }
  // Init Hcore
  cout << "Init and compute Hcore" << endl;
  Hcore hcore(MNDOpara,infoAOs,all2CenterIntegral);
  hcore.CompueHcoreMatrix(hcoreMatrix);

  ScreenUtils::PrintMatrixNxNSymmetric(infoAOs.orbital.size(),hcoreMatrix);

  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  cout << "Dealloc array: hcoreMatrix" << endl;
  Hcore::Dealloc4HcoreMatrix(hcoreMatrix);
  return EXIT_SUCCESS;
}


