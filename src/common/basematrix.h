/*
   Class Base marix for all other class that require a use of big matrix
*/
#ifndef _BASEMATRIX_H_
#define _BASEMATRIX_H_

#include <string>
using std::string;

class BaseMatrix {
 public:
  BaseMatrix();
  BaseMatrix(size_t);
/***************************************************************************************/ 
  // Variables
  size_t matrixSize_;
  // Matrix container
  double* matrixHold_;
/***************************************************************************************/ 
  // Methods
  bool Alloc4Matrix(string matrixName);
  bool Dealloc4Matrix();
  void ComputeMatrix();
/***************************************************************************************/ 
 protected:
  virtual double ComputeElementMatrix(const size_t &i,const size_t &j);
/***************************************************************************************/ 
 private:
  void AssignValue2Matrix(const size_t &i,const size_t &j,const double &value);
};

#endif // _BASEMATRIX_H_ 
