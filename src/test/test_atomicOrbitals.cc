#include <cstdlib>
#include <iomanip>
#include <iostream>
using std::cout;
using std::endl;
using std::setw;
using std::cerr;
#include <string>
using std::string;
#include <vector>
using std::vector;

#include "screenutils.h"
#include "readxyzfile.h"
#include "atom.h"


#include "atomicOrbitals.h"

int main (int argc, char *argv[]){
  cout << endl << "********************************************************" << endl;
  cout << " Testing for  Class AtomicOrbital  " << endl;
  cout << "              Class ListAtomicOrbitals " << endl;
  cout << "********************************************************" << endl << endl;

  string fileName = "filetest_numbers_xyz.xyz";

  cout << "File for read: "<<fileName<<endl;

  ReadXYZFile reader;
  vector<Atom> molecule;

  bool statusAllData = reader.GetValuesFromFile(fileName.c_str(),molecule);

  if (statusAllData) {
    cout << "List of atoms " << endl;
    for (Atom atom : molecule){
      cout << atom.atomSymbol << " " ;
    }
    cout << endl ;
    ListAtomicOrbitals infoAOs;
    infoAOs.SetOrbitals(molecule);

    cout << "Size of infoAOs " << infoAOs.orbital.size() << endl;

    cout << "indexAO  indexAtom  element  angular momentum  int     Coordinates" << endl; 
    for (unsigned int i=0; i< infoAOs.orbital.size() ; i++) {
      cout << setw(4) << infoAOs.orbital[i].indexAO ;
      cout << setw(9) << infoAOs.orbital[i].indexAtom ; 
      cout << setw(12) << infoAOs.orbital[i].element ; 
      cout << setw(9) << "   { " << infoAOs.orbital[i].angularMomentum[0];
      cout << ", " << infoAOs.orbital[i].angularMomentum[1] ; 
      cout << ", " << infoAOs.orbital[i].angularMomentum[2] << " }" ;
      cout << setw(6) << infoAOs.orbital[i].angularMomentumInt;
      cout << setw(7) << "   { " << infoAOs.orbital[i].coordinates[0];
      cout << ", " << infoAOs.orbital[i].coordinates[1] ; 
      cout << ", " << infoAOs.orbital[i].coordinates[2] << " }" << endl; 

    }
    return EXIT_SUCCESS;
  } else {
    cout << "Problem to read file" << endl;
    return EXIT_FAILURE;
  }
}




