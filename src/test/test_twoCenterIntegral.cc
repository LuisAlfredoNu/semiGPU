#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <iomanip>
using std::setw;
using std::setprecision;
#include <string>
using std::string;
#include <vector>
using std::vector;


#include "readxyzfile.h"
#include "atom.h"
#include "atomicOrbitals.h"
#include "MNDO_parameters.h"
#include "fileutils.h"
#include "mymath.h"
#include "screenutils.h"

#include "twocenterintegral.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for Class TwoCenterIntegral " << endl;
	cout << "********************************************************" << endl << endl;


  // Make singles orbitals

  AtomicOrbital orbitalA,orbitalB,orbitalC,orbitalD;

  double coorA[3] = {-1.0,-2.0,-3.0};
  int angularMomA[3] = {1,0,0};
  orbitalA.SetElement( 6);
  orbitalA.SetCoordinates(coorA);
  orbitalA.SetAngularMomentum(angularMomA);
  orbitalA.SetIndexAtom(0);

  int angularMomB[3] = {1,0,0};
  orbitalB.SetElement( 6);
  orbitalB.SetCoordinates(coorA);
  orbitalB.SetAngularMomentum(angularMomB);
  orbitalB.SetIndexAtom(0);

  double coorC[3] = {1.0,1.0,1.5};
  int angularMomC[3] = {1,0,0};
  orbitalC.SetElement( 6);
  orbitalC.SetCoordinates(coorC);
  orbitalC.SetAngularMomentum(angularMomC);
  orbitalC.SetIndexAtom(1);

  int angularMomD[3] = {1,0,0};
  orbitalD.SetElement( 6);
  orbitalD.SetCoordinates(coorC);
  orbitalD.SetAngularMomentum(angularMomD);
  orbitalD.SetIndexAtom(1);

  cout << "orbitalA.angularMomentumInt = " << orbitalA.angularMomentumInt << endl;
  cout << "orbitalB.angularMomentumInt = " << orbitalB.angularMomentumInt << endl;
  cout << "orbitalC.angularMomentumInt = " << orbitalC.angularMomentumInt << endl;
  cout << "orbitalD.angularMomentumInt = " << orbitalD.angularMomentumInt << endl;
  cout << endl;
  // Init MNDO parameters
  MNDOparameter MNDOpara;

  // Init TwoCenterIntegral
  TwoCenterIntegral twoCIntegral(MNDOpara);

  double result =  twoCIntegral.ComputeTwoCenterIntegral(orbitalA,orbitalB,orbitalC,orbitalD);
  cout << "result = " << result << endl;
/* */ 

  string moleculeTest = "C2H6";
  vector<Atom> molecule;

  if (moleculeTest == "C2") {
    molecule.push_back(Atom());
    molecule.push_back(Atom());
    double coorAA[3] = { 0.0, 0.0,-3.0};
    double coorCC[3] = { 1.0, 1.0, 1.0};

    molecule[0].setCoordinates(coorAA);
    molecule[0].setAtomNumber(6);
    molecule[1].setCoordinates(coorCC);
    molecule[1].setAtomNumber(6);
  }else if (moleculeTest == "CH4"){
    string fileName = "filetest_ch4.xyz";
      cout << "File for read: "<<fileName<<endl;
      ReadXYZFile reader;
      bool statusAllData = reader.GetValuesFromFile(fileName.c_str(),molecule);
  }else if (moleculeTest == "C2H6"){
    string fileName = "filetest_c2h6.xyz";
      cout << "File for read: "<<fileName<<endl;
      ReadXYZFile reader;
      bool statusAllData = reader.GetValuesFromFile(fileName.c_str(),molecule);
  }

  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);

  double**** all2CenterIntegral = NULL;
  if (TwoCenterIntegral::Alloc4AllTwoCenterIntegral(molecule,all2CenterIntegral)){
    cout << "Correct Alloc" << endl;
  }else{
    cout << "Bad Alloc" << endl;
  }

  twoCIntegral.ComputeAllTwoCenterIntegral(infoAOs,all2CenterIntegral);

  cout << "Get Compute Data" << endl;
  vector<double> computeData;
  for (int i=0;i<molecule.size();i++) {
    for (int j=0;j<=i;j++) {

      if (molecule[i].atomNumber == 1 && molecule[j].atomNumber  == 1) {
        for (int k=0;k<1;k++) {
          for (int l=0;l<1;l++) {
            computeData.push_back(all2CenterIntegral[i][j][k][l]);
          }
        }
        continue;
      }

      if (molecule[i].atomNumber > 1 && molecule[j].atomNumber > 1) {
        for (int k=0;k<10;k++) {
          for (int l=0;l<10;l++) {
            computeData.push_back(all2CenterIntegral[i][j][k][l]);
          }
        }
        continue;
      }
      if (molecule[i].atomNumber > 1 && molecule[j].atomNumber == 1) {
        for (int k=0;k<1;k++) {
          for (int l=0;l<10;l++) {
            computeData.push_back(all2CenterIntegral[i][j][k][l]);
          }
        }
        continue;
      }
    }
  }

  vector<vector<string>> dataCSV;
  bool readFile = FileUtils::ReadCSV("c2_scf_debug_twoInt.csv",dataCSV);

  if (!readFile) {
    cout << "Problem to read CSV" << endl;
    return EXIT_FAILURE;
  }

  vector<double> refData;

  cout << "Get CSV Data" << endl;
  for (auto& row : dataCSV) {
    for (auto& col : row ) {
      refData.push_back(std::stod(col));
    }
  }


  cout << "Compare data" << endl;
  vector<string> listMOPAC = {"SS","SX","XX","SY","XY","YY","SZ","XZ","YZ","ZZ"};

  int decimals=4;
  cout << std::fixed << setprecision(decimals);
  cout << setw(decimals + 8);
  for (auto label : listMOPAC) {
    cout << label << setw(decimals + 16);
  }

  cout << setw(-18) ;
  int j = 0;
  int totalErrors = 0;
	for (int i=0;i<computeData.size();++i) {
    if (i % 10 == 0) {
      cout  << endl;
      cout << listMOPAC[j];
      ++j;
      if (j == 10) {
        j=0;
      }
    }
    if (sameReal(computeData[i],refData[i],1.0e-4)) {
      cout << setw(decimals + 4) << computeData[i] <<setw(2) << "O" << setw(decimals +4) << refData[i];
    }else{
      cout << setw(decimals + 4) << computeData[i];
      ScreenUtils::SetScrRedBoldFont();
      cout<<setw(2) << "X" << setw(decimals +4) << refData[i];
      ScreenUtils::SetScrNormalFont();
      ++totalErrors;
    }
    cout << setw(2)<< "|";
    if ((1+i)%100 == 0) {
      cout << endl << "total Errors = " << totalErrors <<endl;
      totalErrors = 0;
    }
  }
  cout  << endl;
/*  */
 
  cout << "Test for all possible bielectronic  integral" << endl;
  cout << "get 0,0,4,4 = ";
  cout << TwoCenterIntegral::GetValueFromArray(\
      infoAOs.orbital[0],infoAOs.orbital[0],\
      infoAOs.orbital[4],infoAOs.orbital[4],all2CenterIntegral) << endl;

  int nAOs = infoAOs.orbital.size();
  cout << "orbital size = " << nAOs << endl;

  for (int i=0;i<nAOs;i++) {
    for (int j=0;j<nAOs;j++) {
      for (int k=0;k<nAOs;k++) {
        cout << " element_atomindex " << endl;
        for (int l=0;l<nAOs;l++) {
          cout << infoAOs.orbital[i].element << "_" << infoAOs.orbital[i].indexAtom;
          cout << "  " ;
          cout << infoAOs.orbital[j].element << "_" << infoAOs.orbital[j].indexAtom;
          cout << "  " ;
          cout << infoAOs.orbital[k].element << "_" << infoAOs.orbital[k].indexAtom;
          cout << "  " ;
          cout << infoAOs.orbital[l].element << "_" << infoAOs.orbital[l].indexAtom;
          cout << "  " ;
          cout << "[ " << setw(2) <<  i ;
          cout << " "  << setw(2) <<  j ;
          cout << " | "<< setw(2) <<  k ;
          cout << " "  << setw(2) <<  l << setw(3) << " ] = "; 
          cout << TwoCenterIntegral::GetValueFromArray(infoAOs.orbital[i],\
              infoAOs.orbital[j],infoAOs.orbital[k],infoAOs.orbital[l],\
              all2CenterIntegral) << endl; 
        }
      }
    }
  }

  cout << "Dealloc array of all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  return EXIT_SUCCESS;
}


