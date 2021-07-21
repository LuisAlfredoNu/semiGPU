/*
   solmath.cpp

   Created by Juan Manuel Solano Altamirano on 11/03/13.
   e-mail: jmsolanoalt@gmail.com
   This program was developed in The University of Guelph,
   50 Stone Road West, Guelph
   ON, N1G 2W1, Canada
 
 */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "mymath.h"

#define MAXSAMPLESFORBOYSFUNCTION 13

#ifndef DEBUG
#define DEBUG 0
#endif

/* ************************************************************************** */
int factorial(const int n) {
   return ((n == 1 || n == 0) ? 1 : factorial(n - 1) * n);
}
/* ************************************************************************** */
int doubfact(const int n) {
   return ((n == 0 || n == -1) ? 1 : doubfact(n - 2) * n);
}
/* ************************************************************************** */
double magV3(double (&v)[3]) {
   double m=EPSTODIVIDE+(v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]);
   return sqrt(m);
}
/* ************************************************************************** */
void normalizeV3(double (&v)[3]) {
   double m=EPSTODIVIDE+(v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]);
   m=sqrt(m);
   for (int i=0; i<3; i++) {v[i]/=m;}
   return;
}
/* ************************************************************************** */
void crossProductV3(double (&a)[3],double (&b)[3],double (&c)[3]) {
   c[0]=a[1]*b[2]-a[2]*b[1];
   c[1]=a[2]*b[0]-a[0]*b[2];
   c[2]=a[0]*b[1]-a[1]*b[0];
   return;
}
/* ************************************************************************** */
double dotProductV3(double (&a)[3],double (&b)[3]) {
   return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}
/***************************************************************************************/ 
double distancePointsV3(double (&a)[3],double (&b)[3]){
  double r = 0.0e-10;
  for (int i=0;i<3;i++) {
    r += a[i]-b[i];
    r *= r;
  }
  return sqrt(r);
}
/* ************************************************************************** */
double detM3x3(double (&oM)[3][3]) {
   double res=oM[0][0]*(oM[1][1]*oM[2][2]-oM[2][1]*oM[1][2]);
   res+=(oM[0][1]*(oM[2][0]*oM[1][2]-oM[1][0]*oM[2][2]));
   res+=(oM[0][2]*(oM[1][0]*oM[2][1]-oM[2][0]*oM[1][1]));
   return res;
}
/* ************************************************************************** */
void invertM3x3(double (&oM)[3][3],double (&rM)[3][3]) {
   double det=detM3x3(oM);
   if ( det==0.0e0 ) {
      std::cout << "Warning: Det of the matrix == 0.0! Returning zero Matrix."\
                << std::endl;
      for ( int i=0 ; i<3 ; ++i ) {
         for ( int j=0 ; j<3 ; ++j ) {rM[i][j]=0.0e0;}
      }
      return;
   }
   det=1.0e0/det;
   rM[0][0]=-oM[1][2]*oM[2][1]+oM[1][1]*oM[2][2];
   rM[0][1]=oM[0][2]*oM[2][1]-oM[0][1]*oM[2][2];
   rM[0][2]=-oM[0][2]*oM[1][1]+oM[0][1]*oM[1][2];
   rM[1][0]=oM[1][2]*oM[2][0]-oM[1][0]*oM[2][2];
   rM[1][1]=-oM[0][2]*oM[2][0]+oM[0][0]*oM[2][2];
   rM[1][2]=oM[0][2]*oM[1][0]-oM[0][0]*oM[1][2];
   rM[2][0]=-oM[1][1]*oM[2][0]+oM[1][0]*oM[2][1];
   rM[2][1]=oM[0][1]*oM[2][0]-oM[0][0]*oM[2][1];
   rM[2][2]=-oM[0][1]*oM[1][0]+oM[0][0]*oM[1][1];
   for ( int i=0 ; i<3 ; ++i ) {
      for ( int j=0 ; j<3 ; ++j ) { rM[i][j]*=det; }
   }
   return;
}
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */

