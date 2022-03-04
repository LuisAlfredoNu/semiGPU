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
  cout << endl;
  cout << "********************************************************" << endl;
  cout << " Testing for:  Class AtomicOrbital  " << endl;
  cout << "               Class ListAtomicOrbitals " << endl;
  cout << "********************************************************" << endl << endl;

  string fileName = "../Geometries4Test/benzene.xyz";
  
  if (argc > 1) {
    fileName = argv[1];
  }

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

    cout << "Size of infoAOs " << infoAOs.size() << endl;

    cout << "indexAO  indexAtom  element  angular momentum  int  Core   Coordinates" << endl;
    for (unsigned int i=0; i< infoAOs.size() ; i++) {
      cout << setw(4) << infoAOs.orbital[i].indexAO ;
      cout << setw(9) << infoAOs.orbital[i].indexAtom ; 
      cout << setw(12) << infoAOs.orbital[i].element ; 
      cout << setw(9) << "   { " << infoAOs.orbital[i].angularMomentum[0];
      cout << ", " << infoAOs.orbital[i].angularMomentum[1] ; 
      cout << ", " << infoAOs.orbital[i].angularMomentum[2] << " }" ;
      cout << setw(6) << infoAOs.orbital[i].angularMomentumInt;
      cout << setw(6) << infoAOs.orbital[i].GetCoreCharge();
      cout << setw(7) << "   { " << infoAOs.orbital[i].coordinates[0];
      cout << ", " << infoAOs.orbital[i].coordinates[1] ; 
      cout << ", " << infoAOs.orbital[i].coordinates[2] << " }" << endl; 

    }


#ifdef OPENACC_AVIL
    #pragma acc parallel loop present(infoAOs[0:1])
      for (unsigned int i=0; i< infoAOs.size() ; i++) {
        printf ("%4d"   ,infoAOs.orbital[i].indexAO );
        printf ("%9d"   ,infoAOs.orbital[i].indexAtom) ; 
        printf ("%12d"  ,infoAOs.orbital[i].element ); 
        printf ("%9s %d","   { " , infoAOs.orbital[i].angularMomentum[0]);
        printf ("%s %d" , ", ", infoAOs.orbital[i].angularMomentum[1] ); 
        printf ("%s %d" , ", ", infoAOs.orbital[i].angularMomentum[2] );
        printf (" }");
        printf ("%6d"   , infoAOs.orbital[i].angularMomentumInt);
        printf ("%4d"   , infoAOs.orbital[i].GetCoreCharge());
        printf ("%7s %f", "   { " , infoAOs.orbital[i].coordinates[0]);
        printf ("%s %f" , ", ", infoAOs.orbital[i].coordinates[1] ); 
        printf ("%s %f" , ", ", infoAOs.orbital[i].coordinates[2] );
        printf (" }\n");
      }
#endif
    return EXIT_SUCCESS;
  } else {
    cout << "Problem to read file" << endl;
    return EXIT_FAILURE;
  }
}




