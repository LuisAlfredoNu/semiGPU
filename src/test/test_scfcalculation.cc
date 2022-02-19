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
#include "mymath.h"
#include "twocenterintegral.h"
#include "fileutils.h"
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

  if (argc < 3) {
    cout << "Dont input file in arguments " << endl << endl;
    cout << "input files: geometry.xyz MOPAC_data_energies.csv" << endl;
    return EXIT_FAILURE;
  }
/***************************************************************************************/ 
  string filename = argv[1];
  string fileCSV = argv[2];
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

  //cout << "infoAOs :" << endl;
  //for (size_t i=0;i<nAOs;++i) {
  //  cout << "AO "<< setw(3)  << infoAOs.orbital[i].indexAO << " : ";
  //  cout << setw(5) << infoAOs.orbital[i].element;
  //  cout << setw(5) << infoAOs.orbital[i].angularMomentumInt;
  //  cout << setw(10) << infoAOs.orbital[i].coordinates[0];
  //  cout << setw(10) << infoAOs.orbital[i].coordinates[1];
  //  cout << setw(10) << infoAOs.orbital[i].coordinates[2]<< endl;
  //}
  
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

  ScreenUtils::PrintScrStarLine();
  
  SCFprocess.ComputeSCF();

/***************************************************************************************/ 

  ScreenUtils::PrintScrStarLine();
  double finalEnergy = SCFprocess.finalEnergy;

  cout << "Dealloc array: SCFData" << endl;
  SCFprocess.DeallocSCFData();
  cout << "Dealloc array: hcoreMatrix" << endl;
  hcore.Dealloc4Matrix();
  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  

  vector<vector<string>> dataCSV;
  bool readFile = FileUtils::ReadCSV(fileCSV,dataCSV);

  if (!readFile) {
    cout << "Problem to read CSV" << endl;
    return EXIT_FAILURE;
  }
  double refData = std::stod(dataCSV[1][1]);


  ScreenUtils::PrintScrStarLine();
  vector<double> percentError = {0.0001,0.001,0.01,0.1,1.0};

  double limit;
  bool isSame = false;
  for (auto error : percentError) {
    limit = error * refData/ 100.0;
    limit = std::abs(limit);
    if (sameReal(finalEnergy,refData,limit)) {
      isSame = true;
      limit = error;
      break;
    }
  }
  cout << "Compute Electronic Energy = " << finalEnergy << endl;
  cout << "    Ref Electronic Energy = " << refData << endl;

  if (SCFprocess.energyConverge) {
    cout << "                    Error < " << limit << " %" << endl;
  }else{
    cout << "                    Error = " << "NoConverge" << endl;
  }


  ScreenUtils::PrintScrStarLine();

  if (isSame) {
    ScreenUtils::DisplayGreenMessage("Test pass");
    return EXIT_SUCCESS;
  }else{
    ScreenUtils::DisplayErrorMessage("Test fail");
    return EXIT_FAILURE;
  }
}

