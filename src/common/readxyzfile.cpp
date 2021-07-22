/***************************************************************************************/  
/* Class of ReadXYZFile       */
/***************************************************************************************/  
#ifndef _READXYZFILE_CPP_
#define _READXYZFILE_CPP_
#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <cmath>
using std::abs;
/***************************************************************************************/ 
#include "readxyzfile.h"
/***************************************************************************************/  
/***************************************************************************************/  

ReadXYZFile::ReadXYZFile(){ 
	open_without_problems = true;
	Natoms = 0;
	begindata_pos = 0;
}

bool ReadXYZFile::GetValuesFromFile(string filename, vector<Atom> & molecule){

	ifstream XYZFile;

	XYZFile.open(filename.c_str(),ios::in);

	if (!(XYZFile.good())) {
		cout << "Error: File " << filename << " could not be opened...\n";
#if DEBUG
		cout << __FILE__ << ", line: " << __LINE__ << endl;
#endif
		open_without_problems = false;
		return false;
	}

	Natoms = GetNumofAtoms(XYZFile);
   if(!open_without_problems){
      cout << "Error: File " << filename << " does not have the XYZ format..." << endl;
      return false;
   }
	molecule.resize(Natoms);
	GetDataAtoms(XYZFile,molecule);
   //sortingAtoms(molecule);
	open_without_problems = StatusAllData(molecule);

	XYZFile.close();

	return open_without_problems;
}
/***************************************************************************************/ 
int ReadXYZFile::GetNumofAtoms(ifstream &file){

	int numof_Atoms=0;

	file.seekg(0,file.beg);

	file >> numof_Atoms;

	if(numof_Atoms == 0) open_without_problems = false;

	return numof_Atoms;
}
/***************************************************************************************/ 
void ReadXYZFile::GetDataAtoms(ifstream &file, vector<Atom>& molecule){

	int atomnumber=0;
	string symbol;
	vector<double> vectorposition (3);
			
	int i = 0;
	if(TypeDataNumOChar(file)){
		file.seekg(begindata_pos);
		while(i < Natoms ){
			file >> atomnumber;
			file >> vectorposition[0];
			file >> vectorposition[1];
			file >> vectorposition[2];

			molecule[i].setAtomNumber(atomnumber);
			molecule[i].setCoordinates(vectorposition);

			i++;
		}
	}else{
		file.seekg(begindata_pos);
		while(i < Natoms){
			file >> symbol;
			file >> vectorposition[0];
			file >> vectorposition[1];
			file >> vectorposition[2];

			molecule[i].setAtomSymbol(symbol);
			molecule[i].setCoordinates(vectorposition);
			
			i++;
		}
	}
}
/***************************************************************************************/ 
bool ReadXYZFile::TypeDataNumOChar(ifstream &file){

	bool is_number = true;
	int number;
	string line;
	
	file.seekg(0,file.beg);

	for(int i=0; i<2;i++) getline(file,line);
	
	begindata_pos = file.tellg();		

	file >> number;
	if(file.fail()){
		is_number = false;
		file.clear();
	}
	
	return is_number;
}
/***************************************************************************************/ 
bool ReadXYZFile::StatusAllData(vector<Atom> molecule){

	for(int i=0;i<Natoms;i++)
		if(!molecule[i].statusData) return false;

	return true;
}
/***************************************************************************************/ 
void ReadXYZFile::SortingAtoms(vector<Atom>& molecule){

	// Select witch axe sort firts 
	int sizemolecule = molecule.size();
	vector<double> max_r (3,0.0);
	vector<double> min_r (3,0.0);

	for(int i = 0;i < sizemolecule;i++){	
		if(max_r[0] < molecule[i].atomCoordinates[0]) max_r[0] = molecule[i].atomCoordinates[0];
		if(max_r[1] < molecule[i].atomCoordinates[1]) max_r[1] = molecule[i].atomCoordinates[1];
		if(max_r[2] < molecule[i].atomCoordinates[2]) max_r[2] = molecule[i].atomCoordinates[2];
		                                                     
		if(min_r[0] > molecule[i].atomCoordinates[0]) min_r[0] = molecule[i].atomCoordinates[0];
		if(min_r[1] > molecule[i].atomCoordinates[1]) min_r[1] = molecule[i].atomCoordinates[1];
		if(min_r[2] > molecule[i].atomCoordinates[2]) min_r[2] = molecule[i].atomCoordinates[2];
	}

	vector<double> range_r (3,0.0);
	for(int i = 0;i<3;i++) range_r[i] = max_r[i] - min_r[i];

	vector<int> order2compare (3,0);
	double max = range_r[0], min = range_r[0];

	for(int i = 0;i<3;i++){
		if(max < range_r[i]){
			max = range_r[i];
			order2compare[0] = i;
		}
		if(min > range_r[i]){
			min = range_r[i];
			order2compare[2] = i;
		}
	}
	for(int i = 0;i<3;i++){
		if((max >= range_r[i] && i != order2compare[0]) && (range_r[i] >= min && i != order2compare[2]))
			order2compare[1] = i;
	}
/*
	cout << "Order of range " << endl;
	cout << " range 0 = "<< range_r[0] ; 
	cout << " range 1 = "<< range_r[1] ; 
	cout << " range 2 = "<< range_r[2] << endl ;
	cout << " order 0 = "<< order2compare[0];
	cout << " order 1 = "<< order2compare[1];
	cout << " order 2 = "<< order2compare[2] << endl;
*/
	// Here I adapt the Comb sort algorithm from Wikipedia  
	Atom atom4swap_tmp;

	int swaps = 1;
	int gap = molecule.size();
	double epsilonZ = 0.1001;
	double epsilonY = 0.1001;

	while(!(swaps == 0 && gap ==1)){
		if(gap > 1){
			gap /= 1.3;
			if (gap ==10 || gap == 9){
				gap = 11;
			}
		}

		int i=0; swaps=0;
		while(i + gap < sizemolecule){

			//Sorting by Z
			if((molecule[i].atomCoordinates[order2compare[0]] > molecule[i+gap].atomCoordinates[order2compare[0]]) && (abs(molecule[i].atomCoordinates[order2compare[0]] - molecule[i+gap].atomCoordinates[order2compare[0]]) > epsilonZ)){

				atom4swap_tmp = molecule[i];
				molecule[i] = molecule[i+gap];
				molecule[i+gap] = atom4swap_tmp;

				swaps += 1; 
			}else{
				//Sorting by Y
				if(abs(molecule[i].atomCoordinates[order2compare[0]] - molecule[i+gap].atomCoordinates[order2compare[0]]) < epsilonZ ){
					if(molecule[i].atomCoordinates[order2compare[1]] > molecule[i+gap].atomCoordinates[order2compare[1]] && (abs(molecule[i].atomCoordinates[order2compare[1]] - molecule[i+gap].atomCoordinates[order2compare[1]]) > epsilonY)){

						atom4swap_tmp = molecule[i];
						molecule[i] = molecule[i+gap];
						molecule[i+gap] = atom4swap_tmp;

						swaps += 1; 
					}else{
						// Sorting by X
						if(abs(molecule[i].atomCoordinates[order2compare[1]] - molecule[i+gap].atomCoordinates[order2compare[1]]) < epsilonY){
							if(molecule[i].atomCoordinates[order2compare[2]] > molecule[i+gap].atomCoordinates[order2compare[2]] ){

								atom4swap_tmp = molecule[i];
								molecule[i] = molecule[i+gap];
								molecule[i+gap] = atom4swap_tmp;

								swaps += 1; 
							}
						}
					}
				}
			}
			i += 1;
		}
	}
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _READXYZFILE_CPP_

