/*
   Class Base marix for all other class that require a use of big matrix
*/
#ifndef _BASEMATRIX_H_
#define _BASEMATRIX_H_

#include <string>

class BaseMatrix {
 public:
  BaseMatrix();
/***************************************************************************************/ 
  // Variables
/***************************************************************************************/ 
  // Methods
  static bool Alloc4Matrix(string matrixName,size_t Nrows,double* &matrixHold);
  static bool Dealloc4Matrix(double* &matrixHold);
  virtual void ComputeMatrix(double* &);
 protected:
  void AssignValue2Matrix(const size_t &i,const size_t &j,const double &value,\
       double* &matrixHold);
};

#endif // _BASEMATRIX_H_ 
