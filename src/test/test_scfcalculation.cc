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

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for Class SCFcalculation " << endl;
	cout << "********************************************************" << endl << endl;

  // Init MNDO parameters
  MNDOparameter MNDOpara;

  // Init molecule
  string moleculeTest = "C2";

  vector<Atom> molecule;

  if (moleculeTest == "H6") {
    string fileName = "filetest_h6.xyz";
    cout << "File for read: "<<fileName<<endl;
    ReadXYZFile reader;
    bool statusAllData = reader.GetValuesFromFile(fileName.c_str(),molecule);
  }else if (moleculeTest == "C2") {

    molecule.push_back(Atom());
    molecule.push_back(Atom());
    double coorA[3] = { -1.0,-2.0,-3.0};
    double coorC[3] = {1.0,1.0,1.5};
    molecule[0].setCoordinates(coorA[0],coorA[1],coorA[2]);
    molecule[0].setAtomNumber(6);
    molecule[1].setCoordinates(coorC[0],coorC[1],coorC[2]);
    molecule[1].setAtomNumber(6);
  }else if (moleculeTest == "CH4") {
    string fileName = "filetest_ch4.xyz";
    cout << "File for read: "<<fileName<<endl;
    ReadXYZFile reader;
    bool statusAllData = reader.GetValuesFromFile(fileName.c_str(),molecule);
  }else if (moleculeTest == "C2H6") {
    string fileName = "filetest_c2h6.xyz";
    cout << "File for read: "<<fileName<<endl;
    ReadXYZFile reader;
    bool statusAllData = reader.GetValuesFromFile(fileName.c_str(),molecule);
  }

  cout << "Geometry :" << endl;
  for (size_t i=0;i<molecule.size();++i) {
    cout << setw(3)  << molecule[i].atomSymbol ;
    cout << setw(10) << molecule[i].atomCoordinates[0] ;
    cout << setw(10) << molecule[i].atomCoordinates[1] ;
    cout << setw(10) << molecule[i].atomCoordinates[2] << endl ;
  }

  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);
  int nAOs = infoAOs.orbital.size();

  cout << "infoAOs" << endl;
  for (size_t i=0;i<nAOs;++i) {
    cout << "AO "<< setw(3)  << infoAOs.orbital[i].indexAO << " : ";
    cout << setw(5) << infoAOs.orbital[i].element;
    cout << setw(5) << infoAOs.orbital[i].angularMomentumInt;
    cout << setw(10) << infoAOs.orbital[i].coordinates[0];
    cout << setw(10) << infoAOs.orbital[i].coordinates[1];
    cout << setw(10) << infoAOs.orbital[i].coordinates[2]<< endl;
  }
  // Alloc for all 2 center integrals 
  double**** all2CenterIntegral = NULL;
  if (TwoCenterIntegral::Alloc4AllTwoCenterIntegral(molecule,all2CenterIntegral)){
    cout << "Correct Alloc: all2CenterIntegral" << endl;
  }else{
    cout << "Bad Alloc: all2CenterIntegral " << endl;
  }
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
  /*
  hcore.matrixHold_[1] *= -1.0 ;
  hcore.matrixHold_[3] *= -1.0 ;
  hcore.matrixHold_[6] *= -1.0 ;
  */

  cout << "hcore matrix" << endl;
  ScreenUtils::PrintMatrixNxNSymmetric(nAOs,hcore.matrixHold_);

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

  cout << "Dealloc array: SCFData" << endl;
  SCFprocess.DeallocSCFData();
  cout << "Dealloc array: hcoreMatrix" << endl;
  hcore.Dealloc4Matrix();
  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  return EXIT_SUCCESS;
}


