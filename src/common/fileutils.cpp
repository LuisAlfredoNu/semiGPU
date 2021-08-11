#ifndef _FILEUTILS_CPP_
#define _FILEUTILS_CPP_

#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <sstream>
using std::stringstream;

#include "fileutils.h"
/***************************************************************************************/ 
FileUtils::FileUtils() {
}
/***************************************************************************************/ 
bool FileUtils::ReadCSV(const string filename,vector<vector<string>>& data){
  ifstream ifil(filename.c_str());
  if ( !ifil.good() ) {
    cout << "Error while openning the file " << filename << endl;
    ifil.close();
    return false;
  }
  
  vector<string> row;
  string line, s,word;
  
  while (!ifil.eof()) {

    row.clear();
    std::getline(ifil,line);
    if ( line.length() == 0 ) { 
      break;
    }
    // used for breaking words
    stringstream s(line);
    // read every column data of a row and store it in a string variable, 'word'
    while (getline(s, word,',')) {
      // add all the column data of a row to a vector
      row.push_back(word);
    }
    data.push_back(row);
  }
  ifil.close();
  return true;
}

#endif // _FILEUTILS_CPP_
