#include <cstdlib>
#include <iomanip>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "readxyzfile.h"
#include "atom.h"

int main (int argc, char *argv[]){
	cout << endl << "********************************************************" << endl;
	cout << " Testing for ReadXYZFile Class " << endl;
	cout << "********************************************************" << endl << endl;

	if(argc>1){
		cout << "File for read: "<<argv[1]<<endl;

		ReadXYZFile reader;
		vector<Atom> molecule;

		bool statusAllData = reader.GetValuesFromFile(argv[1],molecule);

		string statusanswer;
		for(unsigned int i=0;i<molecule.size();i++){
			cout << "Atom = " << i << endl;
			for (int ii=0;ii<3;ii++) cout << "Coordinate = " << molecule[i].atomCoordinates[ii] << endl;
			cout << "Element: "<< molecule[i].atomSymbol <<"   Atomic number="<< molecule[i].atomNumber << "   Atomic weight= "<< molecule[i].atomWeight<< endl;
			molecule[i].statusData ? statusanswer="yes" : statusanswer="no"; 
			cout << "All data is fine? " << statusanswer << endl;
			cout << "********************************************************" << endl << endl;
		}   
		cout << "All data is fine ?  " << (statusAllData ? statusanswer="yes" : statusanswer="no") << endl;
		for(unsigned int i=0;i<molecule.size();i++){
			cout << "  " << molecule[i].atomNumber;
			for (int ii=0;ii<3;ii++) cout << std::setw(15)<< std::setprecision(6) << molecule[i].atomCoordinates[ii];
			cout << endl;
		}


		return EXIT_SUCCESS;
	}else { 
		cout << "Dont input file in arguments " << endl << endl;
		return EXIT_FAILURE;
	}
}




