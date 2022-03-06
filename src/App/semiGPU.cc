#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::setw;
#include <vector>
using std::vector;
#include <chrono>
/***************************************************************************************/ 
#include "screenutils.h"
#include "readxyzfile.h"
#include "atom.h"
#include "atomicOrbitals.h"
#include "MNDO_parameters.h"
#include "twocenterintegral.h"
#include "hcore.h"
#include "scfcalculation.h"
#include "corecorerepulsion.h"
/***************************************************************************************/ 
int main (int argc, char *argv[]) {

  auto start = std::chrono::high_resolution_clock::now();
  
  string filename_molecule = argv[argc-1];

  cout << "XYZ file : " << filename_molecule << endl;

  // Read xyz file
  ReadXYZFile reader;
  vector<Atom> molecule;
  if ( ! (reader.GetValuesFromFile(filename_molecule,molecule))){
    cout << "Problem to read xyz file " << endl;
    return EXIT_FAILURE;
  }else{
    cout << "Read XYZ file successful" << endl;
  }

  Atom::PrintGeometry(molecule);

  ScreenUtils::PrintScrStarLine();
  
  // Construct all Atomic Orbitals
  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);
  // Start for compute all two center integrals
  // First alloc memory
  double**** all2CenterIntegrals = NULL;
  if( ! (TwoCenterIntegral::Alloc4AllTwoCenterIntegral(molecule,all2CenterIntegrals))){
    cout << "Problem to alloc all2CenterIntegrals" << endl;
    return EXIT_FAILURE;
  }else{
    cout << "Correct Alloc: all2CenterIntegrals" << endl;
  }
  // Init MNDO parameters 
  MNDOparameter parameters;
  // Compute all two center integrals
  TwoCenterIntegral twoCenInt(parameters);
  twoCenInt.ComputeAllTwoCenterIntegral(infoAOs,all2CenterIntegrals);
  TwoCenterIntegral::To_device(molecule,all2CenterIntegrals);

  // Start for compute Hcore
  Hcore hcore(parameters,infoAOs,all2CenterIntegrals);
  // Alloc Hcore Matrix
  if (hcore.Alloc4Matrix("hcoreMatrix")){
    cout << "Correct Alloc: hcoreMatrix" << endl;
  }else{
    cout << "Bad Alloc: hcoreMatrix " << endl;
  }
  hcore.ComputeMatrix();
  // Init SCF process 
  SCFCalculation SCFprocess(infoAOs,all2CenterIntegrals,hcore);

  bool goodAllocSCF = SCFprocess.AllocSCFData();
  if (goodAllocSCF) {
    cout << "Correct Alloc: SCFData" << endl;
  }else{
    cout << "Bad Alloc: SCFData " << endl;
    return EXIT_FAILURE;
  }
  
  ScreenUtils::PrintScrStarLine();

  SCFprocess.ComputeSCF();

  double electronicEnergy = SCFprocess.finalEnergy;
  ScreenUtils::PrintScrStarLine();
  cout << "  Electronic Energy = " << setw(10) << electronicEnergy << " eV" << endl;
  ScreenUtils::PrintScrStarLine();

  double coreCoreRepulsionEnergy = CoreCoreRepulsion::ComputeRepulsion(parameters,infoAOs,all2CenterIntegrals);

  cout << "Core-Core Repulsion = " << setw(10) << coreCoreRepulsionEnergy << " eV" << endl;
  ScreenUtils::PrintScrStarLine();

  cout << "       Total Energy = " << setw(10) << electronicEnergy + coreCoreRepulsionEnergy << " eV" << endl;
  ScreenUtils::PrintScrStarLine();

  cout << "Dealloc array: SCFData" << endl;
  SCFprocess.DeallocSCFData();

  cout << "Dealloc array: hcoreMatrix" << endl;
  hcore.Dealloc4Matrix();
  cout << "Dealloc array: all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegrals);

  auto end = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  ScreenUtils::PrintScrStarLine();
  double duration_sec = static_cast<double> ( duration.count());
  duration_sec /= 1000.0;
  cout << "Time taken for calculation: " << duration_sec << " s" << endl;

  return EXIT_SUCCESS;
}


