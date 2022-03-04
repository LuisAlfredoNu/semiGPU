#ifndef _SCFCALCULATION_CPP_
#define _SCFCALCULATION_CPP_

#include <iostream>
using std::cout;
using std::endl;

#include <chrono>

#include "mymemory.h"
#include "mymath.h"
#include "electronicenergy.h"

#include "screenutils.h"

#include "scfcalculation.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
SCFCalculation::SCFCalculation(const ListAtomicOrbitals &infoAOs,\
      double**** &all2CenterIntegral,const Hcore &hcore){
  infoAOs_ = &infoAOs;
  nAOs_ = infoAOs.size();
  all2CIntegral_ = all2CenterIntegral;
  hcore_ = &hcore;
  printInfo = false;
}
/***************************************************************************************/ 
bool SCFCalculation::AllocSCFData(){
  //Create array for eigenvec
  if ( ! MyMemory::Alloc1DRealArray("eigenVec",nAOs_*nAOs_,eigenVec,0.0) ) {
    return false;
  }
  //Create array for eigenval
  if (! MyMemory::Alloc1DRealArray("eigenVal",nAOs_,eigenVal,0.0)) {
    return false;
  }
  Pmatrix = new DensityMatrix(*infoAOs_,eigenVec,nAOs_);
  if ( ! Pmatrix->Alloc4Matrix("Pmatrix") ){
    return false;
  }
  Fmatrix = new FockMatrix(*infoAOs_,*hcore_,*Pmatrix,all2CIntegral_);
  if ( ! Fmatrix->Alloc4Matrix("Pmatrix") ){
    return false;
  }
  return true;
}
/***************************************************************************************/ 
bool SCFCalculation::DeallocSCFData(){
  //Create array for eigenvec
  if ( ! MyMemory::Dealloc1DRealArray(eigenVec) ) {
    return false;
  }
  //Create array for eigenval
  if (! MyMemory::Dealloc1DRealArray(eigenVal) ) {
    return false;
  }
  if ( ! Pmatrix->Dealloc4Matrix() ){
    return false;
  }
  if ( ! Fmatrix->Dealloc4Matrix() ){
    return false;
  }
  return true;
}
/***************************************************************************************/ 
void SCFCalculation::ComputeSCF(){
  InitialGuess();

  double oldElectronicenergy = 0.0e-10 ;
  double electronicEnergy = ElectronicEnergy::ComputeEnergy(nAOs_,*hcore_,*Fmatrix,*Pmatrix);

  cout << "SCF step : " << 0 << "  Energy = " << electronicEnergy << endl;
  
  unsigned int maxSCFSteps = 100;
  unsigned int SCFSteps = 0;

  double thresholdEnergy = 5.0e-3;
  energyConverge = false;

  int statusLAPACK;

  while (! energyConverge && maxSCFSteps > SCFSteps ) {
    auto start = std::chrono::high_resolution_clock::now();

    auto startPartial = std::chrono::high_resolution_clock::now();
    statusLAPACK = GetEigenValVecOfFmatrix();
    if (statusLAPACK > 0) {
      cout << "LAPACK: The algorithm failed to converge" << endl;
      break;
    }else if ( statusLAPACK < 0) {
      cout << "LAPACL: the argument "<< statusLAPACK << " had an illegal value" << endl;
      break;
    }
    auto stopLapack = std::chrono::high_resolution_clock::now();
    auto durationLapack = std::chrono::duration_cast<std::chrono::milliseconds>(stopLapack - startPartial);
    
    startPartial = std::chrono::high_resolution_clock::now();
    Pmatrix->ComputeMatrix();
    auto stopPmatrix = std::chrono::high_resolution_clock::now();
    auto durationPmat = std::chrono::duration_cast<std::chrono::milliseconds>(stopPmatrix - startPartial);
    
    startPartial = std::chrono::high_resolution_clock::now();
    Fmatrix->ComputeMatrix();
    auto stopFmatrix = std::chrono::high_resolution_clock::now();
    auto durationFmat = std::chrono::duration_cast<std::chrono::milliseconds>(stopFmatrix - startPartial);
    
    startPartial = std::chrono::high_resolution_clock::now();
    oldElectronicenergy = electronicEnergy;
    electronicEnergy = ElectronicEnergy::ComputeEnergy(nAOs_,*hcore_,*Fmatrix,*Pmatrix);
    auto stopEnergy = std::chrono::high_resolution_clock::now();
    auto durationEnergy = std::chrono::duration_cast<std::chrono::milliseconds>(stopEnergy - startPartial);

    energyConverge = sameReal(oldElectronicenergy,electronicEnergy,thresholdEnergy);

    if (printInfo) {
      cout << "SCFSteps = " << SCFSteps << "  Eigen Values" << endl;
      ScreenUtils::PrintVectorN(nAOs_,eigenVal);
      cout << "Eigen Vectors" << endl;
      ScreenUtils::PrintMatrixNxN(nAOs_,eigenVec);

      cout << endl;
      cout << "SCF : Pmatrix" << endl;
      ScreenUtils::PrintMatrixNxNSymmetric(nAOs_,Pmatrix->matrixHold_);
      cout << "SCF : Fmatrix" << endl;
      ScreenUtils::PrintMatrixNxNSymmetric(nAOs_,Fmatrix->matrixHold_);

      cout << endl << "Electronic Energy  = " << electronicEnergy << endl;

      ScreenUtils::PrintScrStarLine();
    }
    SCFSteps++;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    cout << "SCF step : " << SCFSteps ;
    cout << "  Energy = " << electronicEnergy;
    cout << "  LapackTime = " << durationLapack.count();
    cout << "  PmatTime = " << durationPmat.count();
    cout << "  FmatTime = " << durationFmat.count();
    cout << "  EnerTime = " << durationEnergy.count();
    cout << "  StepTime = " << duration.count() << endl;
  }
  if (printInfo) {
    cout << "Electronic Energy end  = " << electronicEnergy << endl;
  }
  finalEnergy = electronicEnergy;
}
/***************************************************************************************/ 
void SCFCalculation::InitialGuess(){
  Pmatrix->GuessDensityMatrixSwitch(true);
  Pmatrix->ComputeMatrix();
  Pmatrix->GuessDensityMatrixSwitch(false);

  Fmatrix->ComputeMatrix();
  //for (size_t i=0;i<nAOs_;++i) {
  //  Fmatrix->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,i)] -= 0.5;
  //}
  if (printInfo) {
    cout << "Hcore matrix" << endl;
    ScreenUtils::PrintMatrixNxNSymmetric(nAOs_,hcore_->matrixHold_);
    cout << "Guess : Pmatrix" << endl;
    ScreenUtils::PrintMatrixNxNSymmetric(nAOs_,Pmatrix->matrixHold_);
    cout << "Guess : Fmatrix" << endl;
    ScreenUtils::PrintMatrixNxNSymmetric(nAOs_,Fmatrix->matrixHold_);
  }
}
/***************************************************************************************/ 
// LAPACK funtion
// Ref: http://www.netlib.org/lapack/explore-html-3.4.2/d2/d8a/group__double_s_yeigen.html 
extern "C" void  dsyev_(char* JOBZ,char* UPLO,int* N,double* A,int* LDA,double* W,double* WORK,int* LWORK,int* INFO);
/***************************************************************************************/ 
int SCFCalculation::GetEigenValVecOfFmatrix(){

  //Construct eigenVec from Fmatrix_
  //eigenVec column-major-order
  for (size_t i=0;i<nAOs_;++i) {
    for (size_t j=0;j<=i;++j) {
      eigenVec[j*nAOs_ + i] = Fmatrix->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)];
    }
  }

  // Vars 4 lapack
  int matrixOrder = nAOs_, lda = matrixOrder;
  char job = 'V';
  char UPLO = 'L';
  int statusLAPACK;
  int lwork = -1;
  double wkopt;
  /* Query and allocate the optimal workspace */
  dsyev_( &job, &UPLO, &matrixOrder, eigenVec, &lda, eigenVal, &wkopt, &lwork, &statusLAPACK );
 
  // Alloc optimal workspace
  lwork = (int)wkopt;
  double* work;
  MyMemory::Alloc1DRealArray("lapack work",lwork,work);

  /* Solve eigenproblem */
  dsyev_( &job, &UPLO, &matrixOrder, eigenVec, &lda, eigenVal, work, &lwork, &statusLAPACK );
  
  // Dealloc optimal workspace
  MyMemory::Dealloc1DRealArray(work);

  return statusLAPACK;
}
/***************************************************************************************/ 
#endif // _SCFCALCULATION_CPP_
