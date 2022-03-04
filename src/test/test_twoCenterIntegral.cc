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

static double GetValueFromCsvMOPAC(const vector<Atom> &molecule,const ListAtomicOrbitals &infoAOs,int pi,int pj,int pk,int pl,vector<double> &dataCSV);

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
  
  ScreenUtils::PrintScrStarLine();
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
/***************************************************************************************/ 

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
  try{
    for (auto& row : dataCSV) {
      for (auto& col : row ) {
        refData.push_back(std::stod(col));
      }
    }
  }catch(...){
  }
/***************************************************************************************/ 
  ScreenUtils::PrintScrStarLine();

  bool findError = false;
  string errorMen = "";

  cout << "Compare data" << endl;
  vector<string> listMOPAC = {"SS","SX","XX","SY","XY","YY","SZ","XZ","YZ","ZZ"};

  int decimals=4;
  cout << std::fixed << setprecision(decimals);
  cout << setw(decimals + 8);
  int totalErrors = 0;
  int globalTotalErrors = 0;
  int countData = 0;
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
        globalTotalErrors += totalErrors;
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
        globalTotalErrors += totalErrors;
        continue;
      }
      if ( (molecule[i].atomNumber > 1 && molecule[j].atomNumber == 1) || (molecule[i].atomNumber == 1 && molecule[j].atomNumber > 1)) {
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
        globalTotalErrors += totalErrors;
        continue;
      }
    }
  }
  if (globalTotalErrors > 0) {
    findError = true;
    errorMen = "Compute the integral\t";
  }
/***************************************************************************************/ 
  ScreenUtils::PrintScrStarLine();
  cout << "Test for all possible bielectronic  integral" << endl;
  
  size_t nAOs = infoAOs.size();
  cout << "orbital size = " << nAOs << endl << endl;

  double refValue,arrayValue;

  globalTotalErrors = 0;
 
  for (size_t pi=0;pi<nAOs;++pi) {
    for (size_t pj=0;pj<nAOs;++pj) {
      for (size_t pk=0;pk<nAOs;++pk) {
        for (size_t pl=0;pl<nAOs;++pl) {
          refValue = GetValueFromCsvMOPAC(molecule,infoAOs,pi,pj,pk,pl,refData);
          arrayValue = TwoCenterIntegral::GetValueFromArray(\
                       infoAOs.orbital[pi],\
                       infoAOs.orbital[pj],\
                       infoAOs.orbital[pk],\
                       infoAOs.orbital[pl],\
                       all2CenterIntegral);
          if ( ! sameReal(refValue,arrayValue,1.0e-4) ) {
            cout << "[" << pi << "][" << pj << "][" << pk << "][" << pl << "] = " ; 

            ScreenUtils::SetScrRedBoldFont();
            cout << setw(decimals + 4) << arrayValue;
            cout<<setw(2) << "X" ;
            ScreenUtils::SetScrNormalFont();
            cout << setw(decimals +4) << refValue << endl;
            globalTotalErrors += 1;
          }
        }
      }
    }
  }
  if (globalTotalErrors > 0) {
    findError = true;
    errorMen += "Retrun integral Value\t";
  }
#ifdef OPENACC_AVIL
  cout << "all2CenterIntegral ptr = " << all2CenterIntegral << endl;
  TwoCenterIntegral::To_device(molecule,all2CenterIntegral);
  #pragma acc parallel loop gang present(infoAOs,all2CenterIntegral)
  for (size_t pi=0;pi<nAOs;++pi) {
    #pragma acc loop worker
    for (size_t pj=0;pj<nAOs;++pj) {
      #pragma acc loop vector
      for (size_t pk=0;pk<nAOs;++pk) {
        #pragma acc loop seq
        for (size_t pl=0;pl<nAOs;++pl) {
          arrayValue = TwoCenterIntegral::GetValueFromArray(\
                       infoAOs.orbital[pi],\
                       infoAOs.orbital[pj],\
                       infoAOs.orbital[pk],\
                       infoAOs.orbital[pl],\
                       all2CenterIntegral);
        }
      }
    }
  }
  if (globalTotalErrors > 0) {
    findError = true;
    errorMen += "Retrun integral Value\t";
  }
#endif
  cout << "Dealloc array of all2CenterIntegral" << endl;
  TwoCenterIntegral::Dealloc4AllTwoCenterIntegral(molecule,all2CenterIntegral);
  
  ScreenUtils::PrintScrStarLine();
  if (findError) {
    ScreenUtils::DisplayErrorMessage(errorMen);
    return EXIT_FAILURE;
  }else{
    ScreenUtils::DisplayGreenMessage("Test Pass");
    return EXIT_SUCCESS;
  }
}
/***************************************************************************************/ 
/***************************************************************************************/ 
static double GetValueFromCsvMOPAC(const vector<Atom> &molecule,const ListAtomicOrbitals &infoAOs,int pi,int pj,int pk,int pl,vector<double> &dataCSV){
  if (infoAOs.orbital[pi].indexAtom != infoAOs.orbital[pj].indexAtom) {
    return 0.0e-16;
  }
  if (infoAOs.orbital[pk].indexAtom != infoAOs.orbital[pl].indexAtom) {
    return 0.0e-16;
  }
  if (infoAOs.orbital[pi].indexAtom < infoAOs.orbital[pk].indexAtom) {
    swapValues(pi,pk);
    swapValues(pj,pl);
  }

  unsigned int countDataCSV = 0;
  for (int i=0;i < infoAOs.orbital[pi].indexAtom ;i++) {
    for (int j=0;j<=i;j++) {

      if (molecule[i].atomNumber == 1 && molecule[j].atomNumber  == 1) {
        countDataCSV += 1;
        continue;
      }
      if (molecule[i].atomNumber > 1 && molecule[j].atomNumber > 1) {
        countDataCSV += 100;
        continue;
      }
      if (molecule[i].atomNumber == 1 && molecule[j].atomNumber > 1) {
        countDataCSV += 10;
        continue;
      }
      if (molecule[i].atomNumber > 1 && molecule[j].atomNumber == 1) {
        countDataCSV += 10;
        continue;
      }
    }
  }
  int lastAtom = infoAOs.orbital[pi].indexAtom;
  for (int k=0;k < infoAOs.orbital[pk].indexAtom ;k++) {

    if (molecule[lastAtom].atomNumber == 1 && molecule[k].atomNumber  == 1) {
      countDataCSV += 1;
      continue;
    }
    if (molecule[lastAtom].atomNumber > 1 && molecule[k].atomNumber > 1) {
      countDataCSV += 100;
      continue;
    }
    if (molecule[lastAtom].atomNumber == 1 && molecule[k].atomNumber > 1) {
      countDataCSV += 10;
      continue;
    }
    if (molecule[lastAtom].atomNumber > 1 && molecule[k].atomNumber == 1) {
      countDataCSV += 10;
      continue;
    }
  }
  int MOPACindex[4][4];
  
  MOPACindex[0][0] = 0;
  MOPACindex[0][1] = 1 ;
  MOPACindex[1][1] = 2 ;
  MOPACindex[0][2] = 3 ;
  MOPACindex[1][2] = 4 ;
  MOPACindex[2][2] = 5 ;
  MOPACindex[0][3] = 6 ;
  MOPACindex[1][3] = 7 ;
  MOPACindex[2][3] = 8 ;
  MOPACindex[3][3] = 9 ;
                     
  if (pi > pj) {
    swapValues(pi,pj);
  }
  if (pk > pl) {
    swapValues(pk,pl);
  }
  
  countDataCSV += MOPACindex[infoAOs.orbital[pi].angularMomentumInt]\
                            [infoAOs.orbital[pj].angularMomentumInt] * 10;
  countDataCSV += MOPACindex[infoAOs.orbital[pk].angularMomentumInt]\
                            [infoAOs.orbital[pl].angularMomentumInt] ;
 
  return dataCSV[countDataCSV];
}
