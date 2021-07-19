/* 
 * Class for Atoms inside of molecules for align and compare 2 molecules
 */
#ifndef _ATOM_H_
#define _ATOM_H_

#include <string>
using std::string;
#include <vector>
using std::vector;
/***************************************************************************************/ 
class Atom{
/***************************************************************************************/ 
public:
	Atom();
/***************************************************************************************/ 
	// Variables
	int atomNumber;
	string atomSymbol;
	vector<double> atomCoordinates; 
	double atomWeight;
  int atomValenceElectrons;
	bool statusData;

/***************************************************************************************/ 
/* Assing values to coordinates X,Y,Z */
	void setCoordinates(double,double,double);
	void setCoordinates(vector<double>);
/* Get the Value of the atom */
	double getXCoordinate();
	double getYCoordinate();
	double getZCoordinate();
/* Assign name and letter to type of element */
	void setAtomSymbol(string);
	void setAtomNumber(int);


/***************************************************************************************/  
/***************************************************************************************/  

protected:
	double xPosition;
	double yPosition;
	double zPosition;

/***************************************************************************************/  
/***************************************************************************************/  

private:
	string convertAtomNumber2AtomSymbol(int);
	int convertAtomSymbol2AtomNumber(string);
/* Assign atomic weight to the element  */
	double setAtomWeight(int);
/* Assign number of electron of valence */
  int setAtomValenceElectrons(int);
};
#endif // _ATOM_H_
