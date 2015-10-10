#ifndef RYDBERGATOM_H
#define RYDBERGATOM_H

#include <vector>
#include <string>

class RydbergAtom
{
    // Header for the RydbergAtom.cpp file
    public:
        RydbergAtom(std::string f);
        virtual ~RydbergAtom();
        double polarizability; // core polarizability
        void printQD(); // Prints the quantum defects
        double getQuantumDefect(int n, int l, double j);
        std::vector<std::vector<double> > quantumDefects; // 2-D vector that will contain the quantum defects
    protected:
    private:
        std::string qdFile;
};

#endif // RYDBERGATOM_H
