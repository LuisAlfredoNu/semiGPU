#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "screenutils.h"
#include "readxyzfile.h"
#include "atom.h"
#include "atomicOrbitals.h"
#include "MNDO_parameters.h"
#include "twocenterintegral.h"

int main (int argc, char *argv[]) {

  string filename_molecule = argv[argc-1];

  // Read xyz file
  ReadXYZFile reader;
  vector<Atom> molecule;
  if ( ! (reader.GetValuesFromFile(filename_molecule,molecule))){
    cout << "Problem to read xyz file " << endl;
    return EXIT_FAILURE;
  }else{
    cout << "Read XYZ file successful" << endl;
  }
  // Construct all Atomic Orbitals
  ListAtomicOrbitals AOs;
  AOs.SetOrbitals(molecule);
  // Start for compute all two center integrals
  // First alloc memory
  double**** all2CenterIntegrals;
  if( ! (TwoCenterIntegral::Alloc4AllTwoCenterIntegral(molecule,all2CenterIntegrals))){
    cout << "Problem to alloc all2CenterIntegrals" << endl;
    return EXIT_FAILURE;
  }else{
    cout << "All memory for two center integral successful" << endl;
  }
  // Init MNDO parameters 
  MNDOparameter parameters;
  // Compute all two center integrals
  TwoCenterIntegral twoCenInt(parameters);
  twoCenInt.ComputeAllTwoCenterIntegral(AOs.orbital,all2CenterIntegrals);

  return EXIT_SUCCESS;
}


