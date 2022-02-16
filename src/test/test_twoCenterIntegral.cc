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

int main (int argc, char *argv[]){
  cout << endl;
  cout << "********************************************************" << endl;
  cout << " Testing for Class TwoCenterIntegral " << endl;
  cout << "********************************************************" << endl << endl;

/***************************************************************************************/ 
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
  TwoCenterIntegral singleTwoCIntegral(MNDOpara);

  double result =  singleTwoCIntegral.ComputeTwoCenterIntegral(orbitalA,orbitalB,orbitalC,orbitalD);
  cout << "result = " << result << endl;

/***************************************************************************************/   
  ScreenUtils::PrintScrStarLine();
  if (argc < 3) {
    cout << "Dont input file in arguments " << endl << endl;
    cout << "input files: geometry.xyz MOPAC_data_2centerIntegral.csv" << endl;
    return EXIT_FAILURE;
  }
  
  string filename = argv[1];
  string  fileCSV = argv[2];
  
  cout << "File for read: " << filename << endl;

  // Init molecule
  ReadXYZFile reader;
  vector<Atom> molecule;

  bool statusAllData = reader.GetValuesFromFile(filename,molecule);

  if (! statusAllData) {
    cout << "Somethig is wrong with xyz file" << endl;
    return EXIT_FAILURE;
  }

  Atom::PrintGeometry(molecule);

/***************************************************************************************/ 

  ListAtomicOrbitals infoAOs;
  infoAOs.SetOrbitals(molecule);

  // Init TwoCenterIntegral
  TwoCenterIntegral twoCIntegral(MNDOpara);

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
      if (molecule[i].atomNumber == 1 && molecule[j].atomNumber > 1) {
        for (int k=0;k<1;k++) {
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
  bool readFile = FileUtils::ReadCSV(fileCSV,dataCSV);

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
  int totalErrors = 0;
  /* for (auto label : listMOPAC) {
     cout << label << setw(decimals + 16);
     }

     cout << setw(-18) ;
     int j = 0;
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
     ScreenUtils::SetScrRedBoldFont();
     cout << setw(decimals + 4) << computeData[i];
     cout<<setw(2) << "X" ;
     ScreenUtils::SetScrNormalFont();
     cout << setw(decimals +4) << refData[i];
     ++totalErrors;
     }
     cout << setw(2)<< "|";
     if ((1+i)%100 == 0) {
     cout << endl << "total Errors = " << totalErrors <<endl;
     totalErrors = 0;
     }
     }
     cout  << endl;
   */
  int countData = 0;
  int countCSV = 0;
  for (int i=0;i<molecule.size();i++) {
    for (int j=0;j<=i;j++) {

      cout << "Atom A : "<< molecule[i].atomSymbol << "  index : " << i;
      cout << " / ";
      cout << "Atom B : "<< molecule[j].atomSymbol << "  index : " << j;
      cout  << endl;
      totalErrors = 0;
      if (molecule[i].atomNumber == 1 && molecule[j].atomNumber  == 1) {
        for (int k=0;k<1;k++) {
          for (int l=0;l<1;l++) {
            if (sameReal(computeData[countData],refData[countData],1.0e-4)) {
              cout << setw(decimals + 4) << computeData[countData] <<setw(2) << "O" << setw(decimals +4) << refData[countData];
            }else{
              ScreenUtils::SetScrRedBoldFont();
              cout << setw(decimals + 4) << computeData[countData];
              cout<<setw(2) << "X" ;
              ScreenUtils::SetScrNormalFont();
              cout << setw(decimals +4) << refData[countData];
              ++totalErrors;
            }
            cout << setw(2)<< "|";
            countData++;
          }
        }
        cout << endl;
        cout << endl << "total Errors = " << totalErrors <<endl;
        continue;
      }

      if (molecule[i].atomNumber > 1 && molecule[j].atomNumber > 1) {
        for (int k=0;k<10;k++) {
          for (int l=0;l<10;l++) {
            if (sameReal(computeData[countData],refData[countData],1.0e-4)) {
              cout << setw(decimals + 4) << computeData[countData] <<setw(2) << "O" << setw(decimals +4) << refData[countData];
            }else{
              ScreenUtils::SetScrRedBoldFont();
              cout << setw(decimals + 4) << computeData[countData];
              cout<<setw(2) << "X" ;
              ScreenUtils::SetScrNormalFont();
              cout << setw(decimals +4) << refData[countData];
              ++totalErrors;
            }
            cout << setw(2)<< "|";
            countData++;
          }
          cout << endl;
        }
        cout << endl;
        cout << endl << "total Errors = " << totalErrors <<endl;
        continue;
      }
      if (molecule[i].atomNumber > 1 && molecule[j].atomNumber == 1 || molecule[i].atomNumber == 1 && molecule[j].atomNumber > 1) {
        for (int k=0;k<1;k++) {
          for (int l=0;l<10;l++) {
            if (sameReal(computeData[countData],refData[countData],1.0e-4)) {
              cout << setw(decimals + 4) << computeData[countData] <<setw(2) << "O" << setw(decimals +4) << refData[countData];
            }else{
              ScreenUtils::SetScrRedBoldFont();
              cout << setw(decimals + 4) << computeData[countData];
              cout<<setw(2) << "X" ;
              ScreenUtils::SetScrNormalFont();
              cout << setw(decimals +4) << refData[countData];
              ++totalErrors;
            }
            cout << setw(2)<< "|";
            countData++;
          }
          cout << endl;
        }
        cout << endl << "total Errors = " << totalErrors <<endl;
        continue;
      }
    }
  }

  cout << "Test for all possible bielectronic  integral" << endl;
  cout << "get 0,0,4,4 = ";
  cout << "[0][0][4][4] = " << TwoCenterIntegral::GetValueFromArray(\
      infoAOs.orbital[0],infoAOs.orbital[0],\
      infoAOs.orbital[4],infoAOs.orbital[4],all2CenterIntegral) << endl;
  cout << "[0][0][4][5] = "<< TwoCenterIntegral::GetValueFromArray(\
      infoAOs.orbital[0],infoAOs.orbital[0],\
      infoAOs.orbital[4],infoAOs.orbital[5],all2CenterIntegral) << endl;
  cout << "[4][5][0][0] = "<< TwoCenterIntegral::GetValueFromArray(\
      infoAOs.orbital[4],infoAOs.orbital[5],\
      infoAOs.orbital[0],infoAOs.orbital[0],all2CenterIntegral) << endl;

  int nAOs = infoAOs.orbital.size();
  cout << "orbital size = " << nAOs << endl << endl;

  vector<int> list_CSV = {102,110};
  vector<int> list_AO_a = {4,1};
  vector<int> list_AO_b = {4,1};
  vector<int> list_AO_c = {1,4};
  vector<int> list_AO_d = {1,4};

  double value_CSV,value_fromArray,value_fromCompute;

  for (size_t i=0;i<list_CSV.size();++i) {
    cout << " [" << list_AO_a[i] ;
    cout << "][" << list_AO_b[i] ;
    cout << "][" << list_AO_c[i] ;
    cout << "][" << list_AO_d[i] << "]" <<endl;
    value_fromArray = TwoCenterIntegral::GetValueFromArray(\
        infoAOs.orbital[list_AO_a[i]],\
        infoAOs.orbital[list_AO_b[i]],\
        infoAOs.orbital[list_AO_c[i]],\
        infoAOs.orbital[list_AO_d[i]],\
        all2CenterIntegral);
    value_fromCompute = twoCIntegral.ComputeTwoCenterIntegral(\
        infoAOs.orbital[list_AO_a[i]],\
        infoAOs.orbital[list_AO_b[i]],\
        infoAOs.orbital[list_AO_c[i]],\
        infoAOs.orbital[list_AO_d[i]]\
        );
    value_CSV = refData[list_CSV[i]];
    cout << "  getFromArray = " << value_fromArray << endl;
    cout << "getFromCompute = " << value_fromCompute << endl;
    cout << "    getFromCSV = " << value_CSV << endl;
    if ( sameReal(value_CSV,value_fromArray,1.0e-4) ) {
      cout << setw(decimals + 4) << value_CSV <<setw(2) << "O" << setw(decimals +4) << value_fromArray;
    }else{
      ScreenUtils::SetScrRedBoldFont();
      cout << setw(decimals + 4) << value_CSV;
      cout<<setw(2) << "X" ;
      ScreenUtils::SetScrNormalFont();
      cout << setw(decimals +4) << value_fromArray;
    }
    cout << setw(2)<< "|" << endl;
  }

  /*
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
     }*/

  cout << "Dealloc array of all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  return EXIT_SUCCESS;
}


