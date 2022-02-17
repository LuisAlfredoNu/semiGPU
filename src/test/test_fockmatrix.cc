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
#include "readxyzfile.h"
#include "atomicOrbitals.h"
#include "MNDO_parameters.h"
#include "mymemory.h"
#include "mymath.h"
#include "twocenterintegral.h"
#include "screenutils.h"
#include "hcore.h"
#include "densitymatrix.h"
#include "fileutils.h"

#include "fockmatrix.h"

int main (int argc, char *argv[]){
  cout << endl;
  cout << "********************************************************" << endl;
  cout << " Testing for Class FockMatrix " << endl;
  cout << "********************************************************" << endl << endl;
  
  if (argc < 3) {
    cout << "Dont input file in arguments " << endl << endl;
    cout << "input files: geometry.xyz MOPAC_data_1Fmatrix.csv" << endl;
    return EXIT_FAILURE;
  }
/***************************************************************************************/ 
  string filename = argv[1];
  string  fileCSV = argv[2];
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
  ScreenUtils::PrintScrStarLine();
/***************************************************************************************/ 
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
  // Init MNDO parameters
  MNDOparameter MNDOpara;
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

  //cout << "hcore matrix" << endl;
  //ScreenUtils::PrintMatrixNxNSymmetric(nAOs,hcore.matrixHold_);

  double* eigenVec;

  if (MyMemory::Alloc1DRealArray("eigenVec",nAOs*nAOs,eigenVec,0.0)) {
    cout << "Correct Alloc: eigenVec" << endl;
  }else{
    cout << "Bad Alloc: eigenVec " << endl;
  }

  for (int i=0;i<nAOs;++i) {
    eigenVec[i * nAOs + i] = 1.0;
  }

  //cout << "Eigen matrix" << endl;
  //ScreenUtils::PrintMatrixNxN(nAOs,eigenVec);


  // Init Pmatrix
  cout << "Init and Compute DensityMatrix" << endl;
  DensityMatrix Pmatrix(eigenVec,nAOs);
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

  //cout << "Pmatrix matrix" << endl;
  //ScreenUtils::PrintMatrixNxNSymmetric(nAOs,Pmatrix.matrixHold_);

  ScreenUtils::PrintScrStarLine();

  cout << "Init and compute FockMatrix" << endl;
  FockMatrix Fmatrix(infoAOs,hcore,Pmatrix,all2CenterIntegral);
  // Alloc Hcore Matrix
  if (Fmatrix.Alloc4Matrix("Fmatrix")){
    cout << "Correct Alloc: Fmatrix" << endl;
  }else{
    cout << "Bad Alloc: Fmatrix " << endl;
  }
  Fmatrix.ComputeMatrix();
  //cout << "Fock Matrix" << endl;
  //ScreenUtils::PrintMatrixNxNSymmetric(nAOs,Fmatrix.matrixHold_); 
/***************************************************************************************/ 
  // Compare date
  ScreenUtils::PrintScrStarLine();
  cout << "Start to compare Fock matrix data" << endl;

  vector<vector<string>> dataCSV;
  bool readFile = FileUtils::ReadCSV(fileCSV,dataCSV);

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
  //ScreenUtils::PrintMatrixNxNSymmetric(infoAOs.orbital.size(),&refData[0]);

  ScreenUtils::PrintScrStarLine();
  int index;
  int totalErrors = 0;
  int decimals=4;
  cout << std::fixed << setprecision(decimals);
  for (size_t k=0;k<infoAOs.orbital.size();k+=4) {
    for (size_t i=1*k;i<infoAOs.orbital.size();++i) {
      for (size_t j=1*k;j<=i;++j) {
        index = MyMemory::GetIndexSymmetricMatrix(i,j);
        if ( sameReal(Fmatrix.matrixHold_[index],refData[index],1.0e-3) ) {
          cout << setw(decimals + 6) << Fmatrix.matrixHold_[index] ;
          cout << setw(2)<< "O" << setw(decimals +6) << refData[index] << " | ";
        }else{
          ScreenUtils::SetScrRedBoldFont();
          cout << setw(decimals + 6) << Fmatrix.matrixHold_[index] ;
          cout <<setw(2) << "X";
          ScreenUtils::SetScrNormalFont();
          cout  << setw(decimals +6) << refData[index];
          cout   << " | ";
          totalErrors++;
        }
        if ((j+1) % 4 == 0) {
          break;
        }
      }
      cout << endl;
    }
    cout  << endl;
    cout  << endl;
  }
  cout << endl << "total Errors = " << totalErrors <<endl;

  ScreenUtils::PrintScrStarLine();

/***************************************************************************************/ 
  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  cout << "Dealloc array: hcoreMatrix" << endl;
  hcore.Dealloc4Matrix();
  cout << "Dealloc array: Pmatrix" << endl;
  Pmatrix.Dealloc4Matrix();
  cout << "Dealloc array: Fmatrix" << endl;
  Fmatrix.Dealloc4Matrix();

  ScreenUtils::PrintScrStarLine();
  if (totalErrors == 0) {
    ScreenUtils::DisplayGreenMessage("Test pass");
    return EXIT_SUCCESS;
  }else{
    ScreenUtils::DisplayErrorMessage("Fail in Fmatrix");
    return EXIT_FAILURE;
  }
}


