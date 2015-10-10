using namespace std;

class StarkHamiltonian
{
public:

    StarkHamiltonian(double stepSize, double rs, double init1, double init2, int minN, int maxN, double MJ, int addedN)
    {
        h = stepSize;
        re = rs;
        x0 = init1;
        x1 = init2;
        nMin = minN;
        nMax = maxN;
        mj = MJ;
        rStop=0;
        addedStates=addedN;
        //cout << h << ", " << re << ", " << x0 << ", " << x1 << "\n";
    }

    double getMatrixElement(double l, double n, double lp, double np, double mj, double jtot, double jtotp, double qd, double qdp, double F);
    int getSizeOfMatrix();
    static const double eFieldConversion = 5.14220632E9;
    static const double energyConversion = -2.0*109737.31864155017;
    int addedStates;

private:
    double h;
    double re;
    double x0;
    double x1;
    double rStop;
    int nMin;
    int nMax;
    double mj;


    double ClebschGordanCoeff(int j1, int j2, int m1, int m2, int j, int m);
    double binom(int n, int r);
    double Wigner3jSymbol(int j1, int j2, int j3, int m1, int m2, int m3);
    int MIN(int n1, int n2);
    int MAX(int n1, int n2);
    int delta(double a, double b);

};
