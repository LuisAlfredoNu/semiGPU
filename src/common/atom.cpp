/* 
 * Class for Atoms inside of molecules for align and compare 2 molecules
 */
#ifndef _ATOM_CPP_
#define _ATOM_CPP_

#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "atom.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

Atom::Atom(){
	atomNumber=0;
	atomSymbol.clear();
	xPosition=0.0;
	yPosition=0.0;
	zPosition=0.0;
	atomWeight=0.0;
  atomValenceElectrons=0;
	statusData=true;
}

void Atom::setCoordinates(double x, double y, double z){
	xPosition=x;
	yPosition=y;
	zPosition=z;
	atomCoordinates[0]=x;
	atomCoordinates[1]=y;
	atomCoordinates[2]=z;
}
void Atom::setCoordinates(vector<double> xyz_coordinates){
	xPosition=xyz_coordinates[0];
	yPosition=xyz_coordinates[1];
	zPosition=xyz_coordinates[2];
	atomCoordinates[0]=xyz_coordinates[0];
	atomCoordinates[1]=xyz_coordinates[1];
	atomCoordinates[2]=xyz_coordinates[3];
}
void Atom::setCoordinates(double (&xyz_coordinates)[3]){
	xPosition=xyz_coordinates[0];
	yPosition=xyz_coordinates[1];
	zPosition=xyz_coordinates[2];
	atomCoordinates[0]=xyz_coordinates[0];
	atomCoordinates[1]=xyz_coordinates[1];
	atomCoordinates[2]=xyz_coordinates[3];
}
/***************************************************************************************/ 
double Atom::getXCoordinate(){
	return xPosition;
}
double Atom::getYCoordinate(){
	return yPosition;
}
double Atom::getZCoordinate(){
	return zPosition;
}
/***************************************************************************************/ 
/***************************************************************************************/ 
void Atom::setAtomSymbol(string element){
	
	atomNumber = convertAtomSymbol2AtomNumber(element);

	if(atomNumber == 0){
		atomSymbol = element;
		atomSymbol +=" Invalid atomic symbol";
		statusData = false;
	}else{
		atomSymbol = element;
		atomWeight = setAtomWeight(atomNumber);
	}
}
void Atom::setAtomNumber(int number){
	
	if(0 < number && number < 109 ){
		atomSymbol = convertAtomNumber2AtomSymbol(number);
		atomNumber = number;
		atomWeight = setAtomWeight(atomNumber);
	}else{
		atomNumber = 0;
		atomSymbol =" Invalid atomic number";
		statusData = false;
	}
}
double Atom::setAtomWeight(int atomnumber){
	double listatomweight[] = {1.008 , 4.003 , 6.941 , 9.012 , 10.811 , 12.011 , 14.007 , 15.999 , 18.998 , 20.180 , 22.990 , 24.305 , 26.982 , 28.086 , 30.974 , 32.065 , 35.453 , 39.948 , 39.098 , 40.078 , 44.956 , 47.867 , 50.942 , 51.996 , 54.938 , 55.845 , 58.933 , 58.693 , 63.546 , 65.390 , 69.723 , 72.640 , 74.922 , 78.960 , 79.904 , 83.800 , 85.468 , 87.620 , 88.906 , 91.224 , 92.906 , 95.940 , 98.000 , 101.070 , 102.906 , 106.420 , 107.868 , 112.411 , 114.818 , 118.710 , 121.760 , 127.600 , 126.905 , 131.293 , 132.906 , 137.327 , 138.906 , 140.116 , 140.908 , 144.240 , 145.000 , 150.360 , 151.964 , 157.250 , 158.925 , 162.500 , 164.930 , 167.259 , 168.934 , 173.040 , 174.967 , 178.490 , 180.948 , 183.840 , 186.207 , 190.230 , 192.217 , 195.078 , 196.967 , 200.590 , 204.383 , 207.200 , 208.980 , 209.000 , 210.000 , 222.000 , 223.000 , 226.000 , 227.000 , 232.038 , 231.036 , 238.029 , 237.000 , 244.000 , 243.000 , 247.000 , 247.000 , 251.000 , 252.000 , 257.000 , 258.000 , 259.000 , 262.000 , 261.000 , 262.000 , 266.000 , 264.000 , 277.000 , 268.000};
	
	return listatomweight[atomnumber-1];
}
int Atom::setAtomValenceElectrons(int n){
  if ( n<=2 ) {
    return n;
  } else if ( n<=10 ) {
    return n-2;
  } else if ( n<=18 ) {
    return n-10;
  } else {
    statusData = false;
  }
  return -1;

}
/***************************************************************************************/  
/***************************************************************************************/  
int Atom::convertAtomSymbol2AtomNumber(string element){
	
	int numelement=0;
	string listatomymbol[] = {"H" ,"He" ,"Li" ,"Be" ,"B" ,"C" ,"N" ,"O" ,"F" ,"Ne" ,"Na" ,"Mg" ,"Al" ,"Si" ,"P" ,"S" ,"Cl" ,"Ar" ,"K" ,"Ca" ,"Sc" ,"Ti" ,"V" ,"Cr" ,"Mn" ,"Fe" ,"Co" ,"Ni" ,"Cu" ,"Zn" ,"Ga" ,"Ge" ,"As" ,"Se" ,"Br" ,"Kr" ,"Rb" ,"Sr" ,"Y" ,"Zr" ,"Nb" ,"Mo" ,"Tc" ,"Ru" ,"Rh" ,"Pd" ,"Ag" ,"Cd" ,"In" ,"Sn" ,"Sb" ,"Te" ,"I" ,"Xe" ,"Cs" ,"Ba" ,"La" ,"Ce" ,"Pr" ,"Nd" ,"Pm" ,"Sm" ,"Eu" ,"Gd" ,"Tb" ,"Dy" ,"Ho" ,"Er" ,"Tm" ,"Yb" ,"Lu" ,"Hf" ,"Ta" ,"W" ,"Re" ,"Os" ,"Ir" ,"Pt" ,"Au" ,"Hg" ,"Tl" ,"Pb" ,"Bi" ,"Po" ,"At" ,"Rn" ,"Fr" ,"Ra" ,"Ac" ,"Th" ,"Pa" ,"U" ,"Np" ,"Pu" ,"Am" ,"Cm" ,"Bk" ,"Cf" ,"Es" ,"Fm" ,"Md" ,"No" ,"Lr" ,"Rf" ,"Db" ,"Sg" ,"Bh" ,"Hs" ,"Mt"};
	
	int i=0;
	bool test=true;
	while (test && i < 109){
		if(element==listatomymbol[i]){
			numelement=i+1;
			test=false;
		}
		i++;
	}
	return numelement;
}
string Atom::convertAtomNumber2AtomSymbol(int number){

	string listatomymbol[] = {"H" ,"He" ,"Li" ,"Be" ,"B" ,"C" ,"N" ,"O" ,"F" ,"Ne" ,"Na" ,"Mg" ,"Al" ,"Si" ,"P" ,"S" ,"Cl" ,"Ar" ,"K" ,"Ca" ,"Sc" ,"Ti" ,"V" ,"Cr" ,"Mn" ,"Fe" ,"Co" ,"Ni" ,"Cu" ,"Zn" ,"Ga" ,"Ge" ,"As" ,"Se" ,"Br" ,"Kr" ,"Rb" ,"Sr" ,"Y" ,"Zr" ,"Nb" ,"Mo" ,"Tc" ,"Ru" ,"Rh" ,"Pd" ,"Ag" ,"Cd" ,"In" ,"Sn" ,"Sb" ,"Te" ,"I" ,"Xe" ,"Cs" ,"Ba" ,"La" ,"Ce" ,"Pr" ,"Nd" ,"Pm" ,"Sm" ,"Eu" ,"Gd" ,"Tb" ,"Dy" ,"Ho" ,"Er" ,"Tm" ,"Yb" ,"Lu" ,"Hf" ,"Ta" ,"W" ,"Re" ,"Os" ,"Ir" ,"Pt" ,"Au" ,"Hg" ,"Tl" ,"Pb" ,"Bi" ,"Po" ,"At" ,"Rn" ,"Fr" ,"Ra" ,"Ac" ,"Th" ,"Pa" ,"U" ,"Np" ,"Pu" ,"Am" ,"Cm" ,"Bk" ,"Cf" ,"Es" ,"Fm" ,"Md" ,"No" ,"Lr" ,"Rf" ,"Db" ,"Sg" ,"Bh" ,"Hs" ,"Mt"};

	return listatomymbol[number-1]; 
}
#endif // _ATOM_CPP_
