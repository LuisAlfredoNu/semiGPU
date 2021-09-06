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
#include "mymath.h"
#include "mymemory.h"
#include "fileutils.h"
#include "screenutils.h"

#include "basematrix.h"
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

  double coorA[3] = { 0.0, 0.0,-3.0};
  double coorC[3] = { 1.0, 1.0, 1.0};
  molecule[0].setCoordinates(coorA[0],coorA[1],coorA[2]);
  molecule[0].setAtomNumber(6);
  molecule[1].setCoordinates(coorC[0],coorC[1],coorC[2]);
  molecule[1].setAtomNumber(6);

  cout << "Coordinates A = " << coorA[0] << "  " << coorA[1] << "  " << coorA[2] << endl;
  cout << "Coordinates B = " << coorC[0] << "  " << coorC[1] << "  " << coorC[2] << endl;

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

  ScreenUtils::PrintMatrixNxNSymmetric(infoAOs.orbital.size(),hcore.matrixHold_);

  vector<vector<string>> dataCSV;
  bool readFile = FileUtils::ReadCSV("c2_scf_debug_Hcore.csv",dataCSV);

  if (!readFile) {
    cout << "Problem to read CSV" << endl;
    return EXIT_FAILURE;
  }

  vector<double> refData;

  cout << "Get CSV Data" << endl;

  for (size_t i=0;i<infoAOs.orbital.size();++i) {
    for (size_t j=0;j<=i;++j) {
      refData.push_back(std::stod(dataCSV[i][j]));
    }
  }
  ScreenUtils::PrintMatrixNxNSymmetric(infoAOs.orbital.size(),&refData[0]);

  int index;
  int totalErrors = 0;
  int decimals=4;
  cout << std::fixed << setprecision(decimals);
  for (size_t i=0;i<infoAOs.orbital.size();++i) {
    for (size_t j=0;j<=i;++j) {
      index = MyMemory::GetIndexSymmetricMatrix(i,j);
      if ( sameReal(hcore.matrixHold_[index],refData[index],1.0e-1) ) {
        cout << setw(decimals + 6) << hcore.matrixHold_[index] ;
        cout << setw(2)<< "O" << setw(decimals +6) << refData[index] << " | ";
      }else{
        ScreenUtils::SetScrRedBoldFont();
        cout << setw(decimals + 6) << hcore.matrixHold_[index] ;
        cout <<setw(2) << "X"; 
        ScreenUtils::SetScrNormalFont();
        cout  << setw(decimals +6) << refData[index];
        cout   << " | ";
        totalErrors++;
      }
    }
    cout << endl;
  }
  cout << endl << "total Errors = " << totalErrors <<endl;

  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  cout << "Dealloc array: hcoreMatrix" << endl;
  hcore.Dealloc4Matrix();
  return EXIT_SUCCESS;
}


