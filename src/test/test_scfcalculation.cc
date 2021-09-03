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
#include "mymemory.h"
#include "twocenterintegral.h"
#include "screenutils.h"
#include "hcore.h"
#include "densitymatrix.h"
#include "fockmatrix.h"

#include "scfcalculation.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for Class SCFcalculation " << endl;
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

  cout << "Coordinates A = " << coorA[0] << "  " << coorA[1] << "  " << coorA[2] << endl;
  cout << "Coordinates B = " << coorC[0] << "  " << coorC[1] << "  " << coorC[2] << endl;

  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);
  int nAOs = infoAOs.orbital.size();

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

  Hcore hcore(MNDOpara,infoAOs,all2CenterIntegral);
  // Alloc Hcore Matrix
  if (hcore.Alloc4Matrix("hcoreMatrix")){
    cout << "Correct Alloc: hcoreMatrix" << endl;
  }else{
    cout << "Bad Alloc: hcoreMatrix " << endl;
  }
  // Init Hcore
  cout << "Init and compute Hcore" << endl;
  hcore.ComputeMatrix();

  cout << "hcore matrix" << endl;
  ScreenUtils::PrintMatrixNxNSymmetric(nAOs,hcore.matrixHold_);
  
  double* firtsGuess;

  if (MyMemory::Alloc1DRealArray("eigenVec",nAOs*nAOs,firtsGuess,0.0)) {
    cout << "Correct Alloc: firtsGuess" << endl;
  }else{
    cout << "Bad Alloc: firtsGuess " << endl;
  }

  for (int i=0;i<nAOs;++i) {
    firtsGuess[i * nAOs + i] = 1.0;
  }

  cout << "firtsGuess matrix" << endl;
  ScreenUtils::PrintMatrixNxN(nAOs,firtsGuess);


  // Init Hcore
  cout << "Init and Compute DensityMatrix" << endl;
  DensityMatrix Pmatrix(firtsGuess,nAOs);
  // Alloc Pmatrix Matrix
  if (Pmatrix.Alloc4Matrix("Pmatrix")){
    cout << "Correct Alloc: Pmatrix" << endl;
  }else{
    cout << "Bad Alloc: Pmatrix " << endl;
  }
  Pmatrix.ComputeMatrix();
  for (int i=0;i<nAOs;++i) {
    Pmatrix.matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,i)] = 1.0;
  }

  cout << "Pmatrix matrix" << endl;
  ScreenUtils::PrintMatrixNxNSymmetric(nAOs,Pmatrix.matrixHold_);

  cout << "Init and compute FockMatrix" << endl;
  FockMatrix Fmatrix(infoAOs,hcore,Pmatrix,all2CenterIntegral);
  // Alloc Hcore Matrix
  if (Fmatrix.Alloc4Matrix("Fmatrix")){
    cout << "Correct Alloc: Fmatrix" << endl;
  }else{
    cout << "Bad Alloc: Fmatrix " << endl;
  }
  Fmatrix.ComputeMatrix();
  cout << "Fock Matrix" << endl;
  ScreenUtils::PrintMatrixNxNSymmetric(nAOs,Fmatrix.matrixHold_); 

  cout << "Init and compute SCF" << endl;
  SCFCalculation SCFcalculation(infoAOs,all2CenterIntegral,hcore,Pmatrix,Fmatrix);

  SCFcalculation.ComputeSCF();

  cout << "EigenVal Matrix" << endl;
  ScreenUtils::PrintVectorN(nAOs,SCFcalculation.eigenVal);
  cout << "EigenVector Matrix" << endl;
  ScreenUtils::PrintMatrixNxN(nAOs,SCFcalculation.eigenVec);

  cout << "Fock Matrix" << endl;
  ScreenUtils::PrintMatrixNxNSymmetric(nAOs,Fmatrix.matrixHold_); 
  
  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  cout << "Dealloc array: hcoreMatrix" << endl;
  hcore.Dealloc4Matrix();
  cout << "Dealloc array: Pmatrix" << endl;
  Pmatrix.Dealloc4Matrix();
  cout << "Dealloc array: Fmatrix" << endl;
  Fmatrix.Dealloc4Matrix();
  return EXIT_SUCCESS;
}


