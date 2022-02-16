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

#include "readxyzfile.h"
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

int main (int argc, char *argv[]){
  cout << endl;
  cout << "********************************************************" << endl;
  cout << " Testing for Class SCFcalculation " << endl;
  cout << "********************************************************" << endl << endl;

  if (argc < 2) {
    cout << "Dont input file in arguments " << endl << endl;
    return EXIT_FAILURE;
  }
/***************************************************************************************/ 
  string filename = argv[1];
  cout << "File for read: " << filename << endl;

  // Init molecule
  ReadXYZFile reader;
  vector<Atom> molecule;

  bool statusAllData = reader.GetValuesFromFile(filename,molecule);

  if (! statusAllData) {
    cout << "Somethig is wrong with xyz file" << endl;
    return EXIT_FAILURE;
  }

  Atom::PrintGeometry(molecule);

/***************************************************************************************/ 

  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);
  int nAOs = infoAOs.orbital.size();

  cout << "infoAOs :" << endl;
  for (size_t i=0;i<nAOs;++i) {
    cout << "AO "<< setw(3)  << infoAOs.orbital[i].indexAO << " : ";
    cout << setw(5) << infoAOs.orbital[i].element;
    cout << setw(5) << infoAOs.orbital[i].angularMomentumInt;
    cout << setw(10) << infoAOs.orbital[i].coordinates[0];
    cout << setw(10) << infoAOs.orbital[i].coordinates[1];
    cout << setw(10) << infoAOs.orbital[i].coordinates[2]<< endl;
  }
  // Alloc for all 2 center integrals 
  ScreenUtils::PrintScrStarLine();
  double**** all2CenterIntegral = NULL;
  if (TwoCenterIntegral::Alloc4AllTwoCenterIntegral(molecule,all2CenterIntegral)){
    cout << "Correct Alloc: all2CenterIntegral" << endl;
  }else{
    cout << "Bad Alloc: all2CenterIntegral " << endl;
  }
  // Init MNDO parameters
  MNDOparameter MNDOpara;
  // Init TwoCenterIntegral
  cout << "Compute all2CenterIntegral" << endl;
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

  ScreenUtils::PrintScrStarLine();
  cout << "Init and compute SCF" << endl;
  SCFCalculation SCFprocess(infoAOs,all2CenterIntegral,hcore);

  bool goodAllocSCF = SCFprocess.AllocSCFData();
  if (goodAllocSCF) {
    cout << "Correct Alloc: SCFData" << endl;
  }else{
    cout << "Bad Alloc: SCFData " << endl;
    return EXIT_FAILURE;
  }

  SCFprocess.ComputeSCF();

  ScreenUtils::PrintScrStarLine();
  cout << "Final Energy = " << SCFprocess.finalEnergy << endl;
  ScreenUtils::PrintScrStarLine();

  cout << "Dealloc array: SCFData" << endl;
  SCFprocess.DeallocSCFData();
  cout << "Dealloc array: hcoreMatrix" << endl;
  hcore.Dealloc4Matrix();
  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  return EXIT_SUCCESS;
}

