#include <cstdlib>
#include <cstring>
#include <cmath>
#include "mymemory.h"

bool MyMemory::Alloc1DRealArray(string ptrname,const int n,double* &thptr) {
   if (!(thptr=new double[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=0.0;}
      return true;
   }
}
bool MyMemory::Alloc1DRealArray(string ptrname,const int n,double* &thptr,const double inval) {
   if (!(thptr=new double[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=inval;}
      return true;
   }
}
bool MyMemory::Alloc1DIntArray(string ptrname,const int n,int* &thptr) {
   if (!(thptr=new int[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=0;}
      return true;
   }
}
bool MyMemory::Alloc1DIntArray(string ptrname,const int n,int* &thptr,const int inval) {
   if (!(thptr=new int[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=inval;}
      return true;
   }
}
bool MyMemory::Alloc1DBoolArray(string ptrname,const int n,bool* &thptr,const bool inval) {
   if (!(thptr=new bool[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DBoolArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=inval;}
      return true;
   }
}
bool MyMemory::Alloc1DBoolArray(string ptrname,const int n,bool* &thptr) {
   return Alloc1DBoolArray(ptrname,n,thptr,false);
}
bool MyMemory::Dealloc1DBoolArray(bool* &tp) {
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
   }
   return true;
}
bool MyMemory::Dealloc1DRealArray(double* &tp) {
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
   }
   return true;
}
bool MyMemory::Dealloc1DIntArray(int* &tp) {
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
   }
   return true;
}
bool MyMemory::Alloc1DStringArray(string ptrname,const int n, string* &thptr) {
   if (!(thptr=new string[n])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc1DStringArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<n; i++) {thptr[i]=' ';}
      return true;
   }
}
bool MyMemory::Dealloc1DStringArray(string* &tp) {
   if (tp!=NULL) {
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }
}
bool MyMemory::Alloc2DRealArray(string ptrname,const int rows,const int cols,double** &thptr) {
   if (!(thptr=new double*[rows])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<rows; i++) {
         if (!(thptr[i]=new double[cols])) {
            std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DRealArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
      }
      for (int i=0; i<rows; i++) {
         for (int j=0; j<cols; j++) {
            thptr[i][j]=0.000000e0;
         }
      }
      return true;
   }
}
bool MyMemory::Alloc2DRealArray(string ptrname,const int rows,const int cols,double** &thptr,const double inval) {
   if (!(thptr=new double*[rows])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<rows; i++) {
         if (!(thptr[i]=new double[cols])) {
            std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DRealArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
      }
      for (int i=0; i<rows; i++) {
         for (int j=0; j<cols; j++) {
            thptr[i][j]=inval;
         }
      }
      return true;
   }
}
bool MyMemory::Dealloc2DRealArray(double** & tp,const int nr) {
   if (tp!=NULL) {
      for (int i=0; i<nr; i++) {
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }
}
bool MyMemory::Alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr) {
   return Alloc2DIntArray(ptrname,rows,cols,thptr,0);
}
bool MyMemory::Alloc2DIntArray(string ptrname,const int rows,const int cols,int** &thptr,const int val) {
   if (!(thptr=new int*[rows])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<rows; i++) {
         if (!(thptr[i]=new int[cols])) {
            std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc2DIntArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
      }
      for (int i=0; i<rows; i++) {
         for (int j=0; j<cols; j++) {
            thptr[i][j]=val;
         }
      }
      return true;
   }
}
bool MyMemory::Dealloc2DIntArray(int** & tp,const int nr) {
   if (tp!=NULL) {
      for (int i=0; i<nr; i++) {
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
   }
   return true;
}
bool MyMemory::Alloc3DRealArray(string ptrname,const int idx1,const int idx2,const int idx3,double*** &thptr) {
   if (!(thptr=new double**[idx1])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc3DRealArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<idx1; i++) {
         if (!(thptr[i]=new double*[idx2])) {
            std::cout << "Warning: cannot allocate "<< ptrname
            <<", in alloc3DRealArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
         else {
            for (int j=0; j<idx2; j++) {
               if (!(thptr[i][j]=new double[idx3])) {
                  std::cout << "Warning: cannot allocate "<< ptrname
                  <<", in Alloc3DRealArray(...) function.\n";
                  std::cout << __FILE__ << "" << __LINE__ << std::endl;
               }
            }
         }
      }
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            for (int k=0; k<idx3; k++) {
               thptr[i][j][k]=0.0e0;
            }
         }
      }
      return true;
   }
}
bool MyMemory::Alloc4DRealArray(string ptrname,const int idx1,const int idx2,\
      const int idx3,const int idx4,double**** &thptr,double val) {
   string errmsg=string("Warning: cannot allocate ")+ptrname\
                 +string(", in Alloc4DRealArray(...) function.\n");
   if (!(thptr=new double***[idx1])) {
      std::cout << errmsg;
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<idx1; i++) {
         if (!(thptr[i]=new double**[idx2])) {
            std::cout << errmsg;
            std::cout << __FILE__ << " " << __LINE__ << std::endl;
         } else {
            for (int j=0; j<idx2; j++) {
               if (!(thptr[i][j]=new double*[idx3])) {
                  std::cout << errmsg;
                  std::cout << __FILE__ << " " << __LINE__ << std::endl;
               } else {
                  for ( int k=0 ; k<idx3 ; k++ ) {
                     if ( !(thptr[i][j][k]=new double[idx4]) ) {
                        std::cout << errmsg;
                        std::cout << __FILE__ << " " << __LINE__ << std::endl;
                     }
                  }
               }
            }
         }
      }
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            for (int k=0; k<idx3; k++) {
               for ( int l=0 ; l<idx4 ; l++ ) {
                  thptr[i][j][k][l]=val;
               }
            }
         }
      }
      return true;
   }
}
bool MyMemory::Dealloc3DRealArray(double*** &tp,const int idx1,const int idx2) {
   if (tp!=NULL) {
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            delete[] tp[i][j];
            tp[i][j]=NULL;
         }
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }
}
bool MyMemory::Dealloc4DRealArray(double**** &tp,const int idx1,const int idx2,\
      const int idx3) {
   if (tp!=NULL) {
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            for ( int k=0 ; k<idx3 ; k++ ) {
               delete[] tp[i][j][k];
               tp[i][j][k]=NULL;
            }
            delete[] tp[i][j];
            tp[i][j]=NULL;
         }
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }
}
bool MyMemory::Alloc3DIntArray(string ptrname,const int idx1,const int idx2,const int idx3,\
      int*** &thptr,const int val) {
   if (!(thptr=new int**[idx1])) {
      std::cout << "Warning: cannot allocate "<< ptrname <<", in Alloc3DIntArray(...) function.\n";
      std::cout << __FILE__ << "" << __LINE__ << std::endl;
      return false;
   } else {
      for (int i=0; i<idx1; i++) {
         if (!(thptr[i]=new int*[idx2])) {
            std::cout << "Warning: cannot allocate "<< ptrname
            <<", in Alloc3DIntArray(...) function.\n";
            std::cout << __FILE__ << "" << __LINE__ << std::endl;
         }
         else {
            for (int j=0; j<idx2; j++) {
               if (!(thptr[i][j]=new int[idx3])) {
                  std::cout << "Warning: cannot allocate "<< ptrname
                  <<", in Alloc3DIntArray(...) function.\n";
                  std::cout << __FILE__ << "" << __LINE__ << std::endl;
               }
            }
         }
      }
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            for (int k=0; k<idx3; k++) {
               thptr[i][j][k]=val;
            }
         }
      }
      return true;
   }
}
bool MyMemory::Dealloc3DIntArray(int*** &tp,const int idx1,const int idx2) {
   if (tp!=NULL) {
      for (int i=0; i<idx1; i++) {
         for (int j=0; j<idx2; j++) {
            delete[] tp[i][j];
            tp[i][j]=NULL;
         }
         delete[] tp[i];
         tp[i]=NULL;
      }
      delete[] tp;
      tp=NULL;
      return true;
   } else {
      return false;
   }
}
bool MyMemory::AppendTo1DRealArray(string ptrname,const int n,double* &thptr,double thenewval) {
   double *tmpptr;
   bool res=Alloc1DRealArray("tmpptr",(n+1),tmpptr);
   if ( res ) {
      //for ( int i=0 ; i<n ; i++ ) {tmpptr[i]=thptr[i];}
      std::memcpy( tmpptr, thptr, n * sizeof(double) );
      tmpptr[n]=thenewval;
      Dealloc1DRealArray(thptr);
      thptr=tmpptr;
   } else {
      std::cout << "Error: something went wrong while trying to append a new"
         << std::endl << "value in array " << ptrname << std::endl;
   }
   return res;
}
bool MyMemory::AllocSymmetricMatrixReal(string ptrname,const int row,\
      double* &thptr,const double val){

  const int totalSize = row * (row+1) / 2 ;
  return Alloc1DRealArray(ptrname,totalSize,thptr,val);
}
int MyMemory::GetIndexSymmetricMatrix(const int row,const int col){
  return row * (row + 1)/ 2 + col;
}
int MyMemory::GetIndexFullSymmetricMatrix(const int row,const int col){
  if (row < col) {
    return col * (col + 1)/ 2 + row;
  }else{
    return row * (row + 1)/ 2 + col;
  }
}
int MyMemory::GetTotalSizeSymmetricMatrix(const int row){
  return row * (row+ 1 ) / 2 ;
}
void MyMemory::GetIndex_ij_SymetricMatrix(const unsigned int index1D,unsigned int (&index2D)[2]){

  float index1Df = static_cast<float>(index1D);
  float rf = 0.5*(-1.0 + sqrtf(1.0 + 8.0 * index1Df));
  unsigned int r = static_cast<unsigned int>(rf) ;
  unsigned int xLowLimit = r*(r+3)/2 +1;
  if (index1D == xLowLimit) {
    r -= 1 ;
  }
  index2D[0] = r;
  index2D[1] = index1D - r*(r-1)/2 - r;
}
bool MyMemory::DeallocSymmetricMatrixReal(double* &thptr){
  return Dealloc1DRealArray(thptr);
}
