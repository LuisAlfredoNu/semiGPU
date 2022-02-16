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
#include "twocenterintegral.h"
#include "overlap.h"
#include "mymath.h"
#include "mymemory.h"
#include "fileutils.h"
#include "screenutils.h"

#include "basematrix.h"
#include "hcore.h"

int main (int argc, char *argv[]){
  cout << endl;
  cout << "********************************************************" << endl;
  cout << " Testing for Class Hcore " << endl;
  cout << "********************************************************" << endl << endl;
  
  if (argc < 3) {
    cout << "Dont input file in arguments " << endl << endl;
    cout << "input files: geometry.xyz MOPAC_data_Hcore.csv" << endl;
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

/***************************************************************************************/ 

  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);

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

  // Init Overlap
  Overlap overlap(infoAOs);
  if (overlap.Alloc4Matrix("overlapMatrix")){
    cout << "Correct Alloc: overlapMatrix" << endl;
  }else{
    cout << "Bad Alloc: overlapMatrix " << endl;
  }
  // Init Hcore
  cout << "Init and compute Overlap" << endl;
  overlap.ComputeMatrix();
  ScreenUtils::PrintMatrixNxNSymmetric(infoAOs.orbital.size(),overlap.matrixHold_);

  ScreenUtils::PrintScrStarLine();

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

/***************************************************************************************/ 
  // Compare data

  ScreenUtils::PrintScrStarLine();
  cout << "Start to compare Hcore matrix data" << endl;

  vector<vector<string>> dataCSV;
  bool readFile = FileUtils::ReadCSV(fileCSV,dataCSV);

  if (!readFile) {
    cout << "Problem to read CSV" << endl;
    return EXIT_FAILURE;
  }

  vector<double> refData;

  ScreenUtils::PrintScrStarLine();
  cout << "Get CSV Data" << endl;

  cout << "data [0][0] = "<<dataCSV[0][0]<<endl;
  for (size_t i=0;i<infoAOs.orbital.size();++i) {
    for (size_t j=0;j<=i;++j) {
      refData.push_back(std::stod(dataCSV[i][j]));
    }
  }
  ScreenUtils::PrintMatrixNxNSymmetric(infoAOs.orbital.size(),&refData[0]);

  ScreenUtils::PrintScrStarLine();
  cout << "Compare Data" << endl;
  int index;
  int totalErrors = 0;
  int decimals=4;
  cout << std::fixed << setprecision(decimals);
  for (size_t k=0;k<infoAOs.orbital.size();k+=4) {
    for (size_t i=1*k;i<infoAOs.orbital.size();++i) {
      for (size_t j=1*k;j<=i;++j) {
        index = MyMemory::GetIndexSymmetricMatrix(i,j);
        if ( sameReal(hcore.matrixHold_[index],refData[index],1.0e-3) ) {
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
  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  cout << "Dealloc array: hcoreMatrix" << endl;
  hcore.Dealloc4Matrix();
  cout << "Dealloc array: overlapMatrix" << endl;
  overlap.Dealloc4Matrix();
  return EXIT_SUCCESS;
}


