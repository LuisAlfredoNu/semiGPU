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

  bool statusAllData = reader.getValuesFromFile(fileName.c_str(),molecule);

  if (statusAllData) {
    cout << "List of atoms " << endl;
    for (Atom atom : molecule){
      cout << atom.atomSymbol << " " ;
    }
    cout << endl ;
    ListAtomicOrbitals infoAOs;
    infoAOs.setOrbitals(molecule);

    cout << "Size of infoAOs " << infoAOs.orbital.size() << endl;

    cout << "indexAO  indexAtom  element  angular momentum" << endl; 
    for (unsigned int i=0; i< infoAOs.orbital.size() ; i++) {
      cout << setw(4) << infoAOs.orbital[i].indexAO ;
      cout << setw(9) << infoAOs.orbital[i].nAtom ; 
      cout << setw(12) << infoAOs.orbital[i].element ; 
      cout << setw(9) << "   { " << infoAOs.orbital[i].angularMomentum[0];
      cout << ", " << infoAOs.orbital[i].angularMomentum[1] ; 
      cout << ", " << infoAOs.orbital[i].angularMomentum[2] << " }" << endl; 

    }
    return EXIT_SUCCESS;
  } else {
    cout << "Problem to read file" << endl;
    return EXIT_FAILURE;
  }
}




