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
#include "mymemory.h"
#include "screenutils.h"

#include "densitymatrix.h"

int main (int argc, char *argv[]){
  cout << endl;
  cout << "********************************************************" << endl;
  cout << " Testing for Class DensityMatrix " << endl;
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
  ScreenUtils::PrintScrStarLine();
/***************************************************************************************/ 
  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);

  double* eigenVec;

  size_t nAOs = infoAOs.size();

  if (MyMemory::Alloc1DRealArray("eigenVec",nAOs*nAOs,eigenVec,0.0)) {
    cout << "Correct Alloc: eigenVec" << endl;
  }else{
    cout << "Bad Alloc: eigenVec " << endl;
  }

  for (size_t i=0;i<nAOs;++i) {
    eigenVec[i * nAOs + i] = 1.0;
  }

  cout << "Eigen matrix" << endl;
  ScreenUtils::PrintMatrixNxN(nAOs,eigenVec);

  DensityMatrix Pmatrix(infoAOs,eigenVec,nAOs);

  // Alloc Hcore Matrix
  if (Pmatrix.Alloc4Matrix("PmatrixMatrix")){
    cout << "Correct Alloc: PmatrixMatrix" << endl;
  }else{
    cout << "Bad Alloc: PmatrixMatrix " << endl;
  }
  cout << "Compute DensityMatrix" << endl;
  Pmatrix.GuessDensityMatrixSwitch(true);
  Pmatrix.ComputeMatrix();
  ScreenUtils::PrintMatrixNxNSymmetric(nAOs,Pmatrix.matrixHold_);
  
  Pmatrix.GuessDensityMatrixSwitch(false);
  Pmatrix.ComputeMatrix();

  ScreenUtils::PrintMatrixNxNSymmetric(nAOs,Pmatrix.matrixHold_);

  cout << "Update eigenVec values" << endl;
  cout << "eigenVec[1]  = -1.0" << endl;
  cout << "eigenVec[14] = -0.5" << endl;
  eigenVec[1] = -1.0;
  eigenVec[14] = -0.5;

  cout << "Compute again DensityMatrix" << endl;
  Pmatrix.ComputeMatrix();

  ScreenUtils::PrintMatrixNxNSymmetric(nAOs,Pmatrix.matrixHold_);
  cout << "Dealloc array: PmatrixMatrix" << endl;
  Pmatrix.Dealloc4Matrix();
  MyMemory::Dealloc1DRealArray(eigenVec);
  return EXIT_SUCCESS;
}


