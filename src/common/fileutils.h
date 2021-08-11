#ifndef _FILEUTILS_H_
#define _FILEUTILS_H_
#include <string>
using std::string;
#include <vector>
using std::vector;

/***************************************************************************************/ 
class FileUtils {
 public:
  FileUtils();
/***************************************************************************************/ 
  static bool ReadCSV(const string filename,vector<vector<string>>& data);
};

#endif /* _FILEUTILS_H_ */
