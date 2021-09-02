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

#include "mymemory.h"
#include "screenutils.h"

#include "densitymatrix.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for Class DensityMatrix " << endl;
	cout << "********************************************************" << endl << endl;

  double* eigenVec;

  int nAOs = 8;

  if (MyMemory::Alloc1DRealArray("eigenVec",nAOs*nAOs,eigenVec,2.0)) {
    cout << "Correct Alloc: eigenVec" << endl;
  }else{
    cout << "Bad Alloc: eigenVec " << endl;
  }

  for (int i=0;i<nAOs;++i) {
    eigenVec[i * nAOs + i] = 2.0;
  }

  cout << "Eigen matrix" << endl;
  ScreenUtils::PrintMatrixNxN(nAOs,eigenVec);

  DensityMatrix Pmatrix(eigenVec,nAOs);

  // Alloc Hcore Matrix
  if (Pmatrix.Alloc4Matrix("PmatrixMatrix")){
    cout << "Correct Alloc: PmatrixMatrix" << endl;
  }else{
    cout << "Bad Alloc: PmatrixMatrix " << endl;
  }
  // Init Hcore
  cout << "Init and compute DensityMatrix" << endl;
  Pmatrix.ComputeMatrix();

  ScreenUtils::PrintMatrixNxNSymmetric(nAOs,Pmatrix.matrixHold_);

  cout << "Update eigenVec values" << endl;
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


