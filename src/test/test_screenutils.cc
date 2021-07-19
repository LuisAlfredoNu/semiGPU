#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include "screenutils.h"

int main (int argc, char *argv[])
{
   ScreenUtils scrut;
   scrut.DisplayErrorMessage("Testing error message.");
   scrut.DisplayWarningMessage("Testing warning message.");
   scrut.SetScrGreenBoldFont();
   cout << "Testing green message." << endl;
   scrut.SetScrNormalFont();
   scrut.PrintScrStarLine();
   scrut.CenterString("Centered message");
   scrut.PrintScrStarLine();
   return EXIT_SUCCESS;
}


