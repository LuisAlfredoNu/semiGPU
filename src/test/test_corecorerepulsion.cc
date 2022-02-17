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

#include "corecorerepulsion.h"

int main (int argc, char *argv[]){
	cout << endl;
  cout << "********************************************************" << endl;
	cout << " Testing for Class CoreCoreRepulsion " << endl;
	cout << "********************************************************" << endl << endl;
/***************************************************************************************/
  if (argc < 3) {
    cout << "Dont input file in arguments " << endl << endl;
    cout << "input files: geometry.xyz MOPAC_data_energies.csv" << endl;
    return EXIT_FAILURE;
  }

  string fileName = argv[1];
  string fileCSV = argv[2];
  cout << "File for read: " << fileName << endl;

  vector<Atom> molecule;
  ReadXYZFile reader;
  bool statusAllData = reader.GetValuesFromFile(fileName.c_str(),molecule);

  if (! statusAllData) {
    cout << "Problems to get geometry info from xyz file" << endl << endl;
    return EXIT_FAILURE;
  }

  Atom::PrintGeometry(molecule);

  // Create all info of Atomic Orbitals
  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);
  /*
  size_t nAOs = infoAOs.orbital.size();

  cout << "infoAOs :" << endl;
  for (size_t i=0;i<nAOs;++i) {
    cout << "AO "<< setw(3)  << infoAOs.orbital[i].indexAO << " : ";
    cout << setw(5)  << infoAOs.orbital[i].element;
    cout << setw(5)  << infoAOs.orbital[i].angularMomentumInt;
    cout << setw(10) << infoAOs.orbital[i].coordinates[0];
    cout << setw(10) << infoAOs.orbital[i].coordinates[1];
    cout << setw(10) << infoAOs.orbital[i].coordinates[2];
    cout << endl;
  }
  */
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
  
  cout << "Compute Core-Core Repulsion" << endl;
  double energyCoreCoreRepulsion = CoreCoreRepulsion::ComputeRepulsion(MNDOpara,infoAOs,all2CenterIntegral);

  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  
  vector<vector<string>> dataCSV;
  bool readFile = FileUtils::ReadCSV(fileCSV,dataCSV);

  if (!readFile) {
    cout << "Problem to read CSV" << endl;
    return EXIT_FAILURE;
  }
  double refData = std::stod(dataCSV[2][1]);


  ScreenUtils::PrintScrStarLine();
  cout << "Compute Core-Core Repulsion = " << energyCoreCoreRepulsion << endl;
  cout << "    Ref Core-Core Repulsion = " << refData << endl;

  ScreenUtils::PrintScrStarLine();

  if (sameReal(energyCoreCoreRepulsion,refData,0.0001)) {
    ScreenUtils::DisplayGreenMessage("Test pass");
    return EXIT_SUCCESS;
  }else{
    ScreenUtils::DisplayErrorMessage("Test fail");
    return EXIT_FAILURE;
  }

}


