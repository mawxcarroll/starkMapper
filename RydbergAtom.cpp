#include "RydbergAtom.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <gsl/gsl_math.h>

using namespace std;

// The RydbergAtom object contains all information about quantum defects and core polarizability. This
// works for any type of atom as long as the quantum defects and core polarizability are entered into a file in the
// appropriate format; see the RubidiumDefects.txt file for an example of this. The aforementioned quantities are necessary
// for the calculation the Stark map.

RydbergAtom::RydbergAtom(string filename)
{
    ifstream ifs(filename.c_str());
    string myLine;
    bool foundPolar = false; // the polarizability should be the first value stored in the file
    stringstream ss; // used to parse lines of defects into vectors
    string defectString; // defects stored in string form
    vector<double> tempDefectsFromFile; // vector of defects for a given state
    // Read until the end of the file
    while(getline(ifs,myLine))
    {
        // ignore comments
        if(myLine[0] != '#')
        {
            // The polarizability should be the first value stored in the file
            if(!foundPolar)
            {
                polarizability = atof(myLine.c_str());
                foundPolar = true;
            }
            else
            {
                // Tokenize the line, convert the strings to doubles, store them in a vector
                ss.str(myLine);
                while(getline(ss,defectString,','))
                    tempDefectsFromFile.push_back(atof(defectString.c_str()));
                // Add the vector for a given state to the total vector
                quantumDefects.push_back(tempDefectsFromFile);
                // clear the temporary vector
                defectString.clear();
                ss.clear();
                tempDefectsFromFile.clear();
            }
        }
    }
    ifs.close();
    //To test, print out each of the defects and the polarizability
    /*cout << "polarizability " << polarizability << endl;
    for(int x=0; x< quantumDefects.size(); x++)
    {
        cout << "defect set number " << x << endl;
        for(int y = 0; y < quantumDefects[x].size(); y++)
        {
            cout << quantumDefects[x][y] << endl;
        }
    }*/
}

RydbergAtom::~RydbergAtom()
{
    //dtor
}

double RydbergAtom::getQuantumDefect(int n, int l, double j)
{
    // The number of defects stored in the file for an atom should be independent of mj,
    // but the smallest possible value of mj is needed for indexing.
    vector<double> defects;
    double result = 0;
    double minMJ = 0.5;
    // Find the index of the state of interest given l and j
    int defectIndex = int(double(l) + j - minMJ);
    // If there is no data for the state, return 0.0
    if(defectIndex >= quantumDefects.size())
    {
        //use the smallest quantum defect to estimate the rest with l^-5 scaling
        //this is basically hard-coded for Rb!
        result = quantumDefects[quantumDefects.size()-1][0]*1024.0/pow((double)l, 5.0);
        return result;
    }
    else
    {
        // If there is data for the state, compute the defect and return the result
        defects = quantumDefects[defectIndex];
        for(int termNumber = 0; termNumber < defects.size(); termNumber++)
        {
            result+= defects[termNumber]/pow(double(n) - defects[0], 2*termNumber);
        }
        return result;
    }
}

void RydbergAtom::printQD()
{
    // Prints all quantum defects corresponding to the same state in one line, then does this until all defects are printed
    for(int i = 0; i < quantumDefects.size();i++)
    {
        for(int j=0; j < quantumDefects[i].size();j++)
        {
            cout << quantumDefects[i][j] << "|||";
        }
        cout << endl;
    }
}

