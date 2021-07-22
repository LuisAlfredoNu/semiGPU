#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "atom.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for Atom Class " << endl;
	cout << "********************************************************" << endl << endl;

	vector<Atom> molecule (4,Atom());

	double x1=1.5, y1=2.0, z1=0.5;
	double x2=3.0, y2=1.0, z2=1.5;
	double x3=1.2, y3=3.4, z3=5.6;
	double x4=2.2, y4=2.4, z4=2.6;
	string nameelemnt01 ("H");
	int numelement02=6;
	string nameelemnt03 ("Alfredo");
	int numelement04=222;

	string statusanswer;

	molecule[0].setCoordinates(x1,y1,z1);
	molecule[0].setAtomSymbol(nameelemnt01);
	molecule[1].setCoordinates(x2,y2,z2);
	molecule[1].setAtomNumber(numelement02);
	molecule[2].setCoordinates(x3,y3,z3);
	molecule[2].setAtomSymbol(nameelemnt03);
	molecule[3].setCoordinates(x4,y4,z4);
	molecule[3].setAtomNumber(numelement04);

	for(unsigned int i=0;i<molecule.size();i++){
		cout << "Atom = " << i << endl;
		for (int ii=0;ii<3;ii++) cout << "Coordinate = " << molecule[i].atomCoordinates[ii] << endl;
		cout << "Element: "<< molecule[i].atomSymbol <<"   Atomic number="<< molecule[i].atomNumber << "   Atomic weight= "<< molecule[i].atomWeight<< endl;
		molecule[i].statusData ? statusanswer="yes" : statusanswer="no"; 
		cout << "All data is fine? " << statusanswer << endl;
		cout << "********************************************************" << endl << endl;
	 }


	return EXIT_SUCCESS;
}


