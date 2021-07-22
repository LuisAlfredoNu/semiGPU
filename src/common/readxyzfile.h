/***************************************************************************************/ 
/*
 * 
 * Class for read a XYZ file  
 * Methods : Public
 * 	getValuesFromFile -> Open/Close the file and pass the ifstream to the other methods  
 * Methods : Private
 * 	getNumofAtoms -> From file read the firts line for kwoning the number of atoms in the molecule
 * 	getDataAtoms -> Assing values from file (atomic number and coordinates) to the object atom
 * 	typeDataNumOChar -> Check if is a number or symbol for atoms
 * 	statusAllData -> Chech if all data are good
 *
 */
/***************************************************************************************/  
#ifndef _READXYZFILE_H_
#define _READXYZFILE_H_
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <string>
using std::string;

#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "atom.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

class ReadXYZFile{
	/***************************************************************************************/ 
	public:
		ReadXYZFile();
	/***************************************************************************************/ 
		// Variables
		int Natoms;
		bool open_without_problems;
	/***************************************************************************************/ 
		bool GetValuesFromFile(string,vector<Atom>&);
		void SortingAtoms(vector<Atom> &);
	
	/***************************************************************************************/
	/***************************************************************************************/ 
	private:
		// Variables
		int begindata_pos;
	/***************************************************************************************/ 
		void GetDataAtoms(ifstream &file, vector<Atom> &);
		int GetNumofAtoms(ifstream &file);
		bool TypeDataNumOChar(ifstream &file);
		bool StatusAllData(vector<Atom>);


	/***************************************************************************************/ 
	/***************************************************************************************/ 

	protected:
};
#endif // _READXYZFILE_H_
