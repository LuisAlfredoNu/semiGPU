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
  BaseMatrix(const BaseMatrix&);
/***************************************************************************************/ 
  // Variables
  // Matrix NxN size
  size_t matrixSize_;
  // Matrix container
  double* matrixHold_;
  // Size of full container
  size_t array1DSize_;
/***************************************************************************************/ 
  // Methods
  bool Alloc4Matrix(string matrixName);
  bool Dealloc4Matrix();
  virtual void ComputeMatrix();
/***************************************************************************************/ 
 protected:
  // Methods
  virtual double ComputeElementMatrix(const size_t &i,const size_t &j);
  void To_deviceMatrix();
  void From_deviceMatrix();
  void Update_hostMatrix();
  void Update_deviceMatrix();
/***************************************************************************************/ 
 private:
  // Methods
  void AssignValue2Matrix(const size_t &i,const size_t &j,const double &value);
};

#endif // _BASEMATRIX_H_ 
