#ifndef _SCFCALCULATION_CPP_
#define _SCFCALCULATION_CPP_

#include "mymemory.h"
#include "electronicenergy.h"

#include <iostream>
using std::cout;
using std::endl;
#include "screenutils.h"

#include "scfcalculation.h"
/***************************************************************************************/ 
/***************************************************************************************/ 
SCFCalculation::SCFCalculation(const ListAtomicOrbitals &infoAOs,\
      double**** &all2CenterIntegral,const Hcore &hcore,\
      DensityMatrix &Pmatrix,FockMatrix &Fmatrix){
  infoAOs_ = &infoAOs;
  nAOs_ = infoAOs.orbital.size();
  all2CIntegral_ = all2CenterIntegral;
  hcore_ = &hcore;
  Pmatrix_ = &Pmatrix;
  Fmatrix_ = &Fmatrix;
  
  //Create array for eigenvec
  goodAllocEigenVec = MyMemory::Alloc1DRealArray("eigenVec",nAOs_*nAOs_,eigenVec,0.0);
  //Create array for eigenval
  goodAllocEigenVal = MyMemory::Alloc1DRealArray("eigenVal",nAOs_,eigenVal,0.0);
}
/***************************************************************************************/ 
void SCFCalculation::ComputeSCF(){
  InitialGuess();
  double electronicEnergy = ElectronicEnergy::ComputeEnergy(nAOs_,*hcore_,*Fmatrix_,*Pmatrix_);
  GetEigenValVecOfFmatrix();
}
/***************************************************************************************/ 
/*
   Plan
   Crear el frirts guess
   Calcular los valores y vectores propios
 */
/***************************************************************************************/ 
void SCFCalculation::InitialGuess(){
  for (size_t i=0;i<nAOs_;++i) {
    Pmatrix_->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,i)] = 1.0;
  }
  Fmatrix_->ComputeMatrix();
  for (size_t i=0;i<nAOs_;++i) {
    Fmatrix_->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,i)] -= 0.5;
  }
}
/***************************************************************************************/ 
extern "C" void  dsyev_(char* JOBZ,char* UPLO,int* N,double* A,int* LDA,double* W,double* WORK,int* LWORK,int* INFO);
/***************************************************************************************/ 
int SCFCalculation::GetEigenValVecOfFmatrix(){

  //Construct eigenVec from Fmatrix_
  //eigenVec column-major-order
  for (size_t i=0;i<nAOs_;++i) {
    for (size_t j=0;j<=i;++j) {
      eigenVec[j*nAOs_ + i] = Fmatrix_->matrixHold_[MyMemory::GetIndexSymmetricMatrix(i,j)];
    }
  }

  // Vars 4 lapack
  int matrixOrder = nAOs_, lda = matrixOrder;
  char job = 'V';
  char UPLO = 'L';
  int info;
  int lwork = -1;
  double wkopt;
  /* Query and allocate the optimal workspace */
  dsyev_( &job, &UPLO, &matrixOrder, eigenVec, &lda, eigenVal, &wkopt, &lwork, &info );
  
  lwork = (int)wkopt;
  double* work;
  MyMemory::Alloc1DRealArray("lapack work",lwork,work);

  /* Solve eigenproblem */
  dsyev_( &job, &UPLO, &matrixOrder, eigenVec, &lda, eigenVal, work, &lwork, &info );
  
  /* Check for convergence */
  if( info > 0 ) {
    printf( "The algorithm failed to compute eigenvalues.\n" );
    exit( 1 );
  }
  MyMemory::Dealloc1DRealArray(work);
}
/***************************************************************************************/ 
#endif // _SCFCALCULATION_CPP_
