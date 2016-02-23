#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include "mpi.h"
#include "netcdf.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "StarkHamiltonian.h"


using namespace std;


double StarkHamiltonian::getMatrixElement(double l, double n, double lp, double np, double mj, double jtot, double jtotp, double qd, double qdp, double F)
{
	//double h; //the step size
	double rs; //the starting point (outer limit)
	//double re; //the ending point (inner limit)
	double x[3], xp[3], g[3], gp[3], r[3], t[3]; //t is the natural log of r
	int i; //for convenience
	//double x0, x1; //the initial conditions
	i = 1; //this is just to make the code readable and similar to the equations
	int lMin = 3;
	double m = 0;
	x[i - 1] = x0;
	x[i] = x1;
	xp[i - 1] = x0;
	xp[i] = x1;
	double element;
	rs = 2 * n * (n + 15.);
	r[i - 1] = rs;
	r[i] = rs * exp(-1.0 * h);
	t[i - 1] = log(r[i - 1]);
	t[i] = log(r[i]);
	double sum1, sum2, sum3;
	double converge[2]; // for convergence testing
	bool converged = false;
	int j = 2;

	double testAnswer = 0;

 //cout << h << ", " << re << ", " << x0 << ", " << x1 << "\n";

	double E = -1 / (2. * (n - qd) * (n - qd));
	double Ep = -1 / (2. * (np - qdp) * (np - qdp));
	//just return the energy if it's a diagonal element
	if (l == lp && n == np && jtot==jtotp)
		return E;
	//calculate the angular part, which is quite simple
	//since there is an exact formula, if it's zero, just return 0
	//need to sum over two values of ml = mj +/- 1/2
	//Also need to multiply by Clebsch-Gordan coefficients
	double angPart = 0.0;
	//printf("hi 2\n");
	//printf("*** j1 =  %d, j2 = %d, m1 = %d, m2 = %d, j = %d, m = %d\n",((int) (l * 2)), 1, (int)(m*2), (int)((mj-m)*2), (int)(jtot*2), (int)(mj*2));
	//printf("hi 2\n");
	//ClebschGordanCalculator cgc = new ClebschGordanCalculator();
	if (lp == (l - 1)) {
		m = mj + .5;
		angPart = sqrt((double) ((l * l - m * m) / ((2 * l + 1) * (2 * l - 1)))) * ClebschGordanCoeff((int) (l * 2), 1, (int) (m * 2), (int) ((mj - m) * 2), (int) (jtot * 2), (int) (mj * 2)) * ClebschGordanCoeff((int) ((l - 1) * 2), 1, (int) (m * 2), (int) ((mj - m) * 2), (int) (jtotp * 2),
				(int) (mj * 2));
		//cout << " xxangpart = "<< angPart << "\n";
		m = mj - .5;
		angPart += sqrt((l * l - m * m) / ((2 * l + 1) * (2 * l - 1))) * ClebschGordanCoeff((int) (l * 2), 1, (int) (m * 2), (int) ((mj - m) * 2), (int) (jtot * 2), (int) (mj * 2)) * ClebschGordanCoeff((int) ((l - 1) * 2), 1, (int) (m * 2), (int) ((mj - m) * 2), (int) (jtotp * 2), (int) (mj * 2));
	}

	if (lp == (l + 1)) {
		m = mj + .5;
		angPart = sqrt(((l + 1) * (l + 1) - m * m) / ((2 * l + 1) * (2 * l + 3))) * ClebschGordanCoeff((int) (l * 2), 1, (int) (m * 2), (int) ((mj - m) * 2), (int) (jtot * 2), (int) (mj * 2)) * ClebschGordanCoeff((int) ((l + 1) * 2), 1, (int) (m * 2), (int) ((mj - m) * 2), (int) (jtotp * 2),
				(int) (mj * 2));
		m = mj - .5;
		angPart += sqrt(((l + 1) * (l + 1) - m * m) / ((2 * l + 1) * (2 * l + 3))) * ClebschGordanCoeff((int) (l * 2), 1, (int) (m * 2), (int) ((mj - m) * 2), (int) (jtot * 2), (int) (mj * 2)) * ClebschGordanCoeff((int) ((l + 1) * 2), 1, (int) (m * 2), (int) ((mj - m) * 2), (int) (jtotp * 2),
				(int) (mj * 2));
	}
	//cout << " angpart = " << angPart << "\n";

	if (angPart == 0)
		return angPart;
	//printf("hi");
	double h2 = pow(h, 2);
	//return 0;
	//calculate all parameters for loop initialization
	r[i + 1] = rs * exp(-2 * h);
	t[i + 1] = log(r[i + 1]);
	g[i - 1] = 2 * exp(2. * t[i - 1]) * (-1 / r[i - 1] - E) + (l + .5) * (l + .5);
	g[i] = 2 * exp(2. * t[i]) * (-1 / r[i] - E) + (l + .5) * (l + .5);
	g[i + 1] = 2 * exp(2. * t[i + 1]) * (-1 / r[i + 1] - E) + (l + .5) * (l + .5);
	gp[i - 1] = 2 * exp(2. * t[i - 1]) * (-1 / r[i - 1] - Ep) + (lp + .5) * (lp + .5);
	gp[i] = 2 * exp(2. * t[i]) * (-1 / r[i] - Ep) + (lp + .5) * (lp + .5);
	gp[i + 1] = 2 * exp(2. * t[i + 1]) * (-1 / r[i + 1] - Ep) + (lp + .5) * (lp + .5);
	sum1 = x[i - 1] * xp[i - 1] * pow(r[i], 3) + x[i] * xp[i] * pow(r[i], 3);
	sum2 = x[i - 1] * r[i - 1] * r[i - 1] + x[i] * r[i] * r[i];
	sum3 = xp[i - 1] * r[i - 1] * r[i - 1] + xp[i] * r[i] * r[i];
	element = converge[i - 1] = sum1 / sqrt(sum2 * sum3);
/*
	 printf("\n");
	 printf("t[i+1] = %f\n", t[i + 1]);
	 printf("E = %f Ep = %f\n", E, Ep);
	 printf("sum1 new = %f\n", sum1);
	 printf("sum2 new = %f\n", sum2);
	 printf("sum3 new = %f\n", sum3);
	 printf("element new = %f\n", element);
	 */
	//return 0;
	if (l > lMin)
		rStop = 0;
	else
		rStop = re;
	//The Numerov Algorithm

	while (r[i + 1] > rStop) {
		//the numerov algorithm calculations
		x[i + 1] = (x[i - 1] * (g[i - 1] - 12. / h2) + x[i] * (10 * g[i] + 24. / h2)) / (12. / h2 - g[i + 1]);
		xp[i + 1] = (xp[i - 1] * (gp[i - 1] - 12. / h2) + xp[i] * (10 * gp[i] + 24. / h2)) / (12. / h2 - gp[i + 1]);
		//calculate contributions to sum
		sum1 += x[i + 1] * xp[i + 1] * pow(r[i + 1], 3);
		sum2 += x[i + 1] * x[i + 1] * r[i + 1] * r[i + 1];
		sum3 += xp[i + 1] * xp[i + 1] * r[i + 1] * r[i + 1];
		//set current points to old points and calculate new points
		x[i - 1] = x[i];
		x[i] = x[i + 1];
		xp[i - 1] = xp[i];
		xp[i] = xp[i + 1];
		g[i - 1] = g[i];
		g[i] = g[i + 1];
		gp[i - 1] = gp[i];
		gp[i] = gp[i + 1];
		j++; //increment loop index
		r[i + 1] = rs * exp(-1 * j * h);
		t[i + 1] = log(r[i + 1]);
		g[i + 1] = 2 * exp(2. * t[i + 1]) * (-1 / r[i + 1] - E) + (l + .5) * (l + .5);
		gp[i + 1] = 2 * exp(2. * t[i + 1]) * (-1 / r[i + 1] - Ep) + (lp + .5) * (lp + .5);
		element = sum1 / sqrt(sum2 * sum3);
		//cout << i << ", " << element << "\n";
		if (!(element >= 0 || element < 0)) {
			element = 0;
			converged = true;
		}

		converge[i - 1] = converge[i];
		converge[i] = element;
		if (converged && l > lMin) {
			if ((fabs(1 - converge[i] / converge[i - 1]) > 1E-6))
				r[i + 1] = rStop;
		}

		if (fabs(1 - converge[i] / converge[i - 1]) < 1E-6 && !converged) {
			converged = true;
			testAnswer = converge[i];
		}
	}

	element *= angPart;
	element *= F;
	return element;
}

//Returns # of rows (and hence, columns) of the Hamiltonian
int StarkHamiltonian::getSizeOfMatrix()
{
    int matrixSize=0; //size of the matrix to be returned
    for(int n=nMin;n<=nMax;n++)
    {
        //for every n state, there are n angular momentum states. These angular momentum states can be spin up or spin down,
        //thus, naively, we say there are 2*n j states. However, for the s state, one of the j states would be j=-1/2,
        //which is non-physical, so we remove this. As a result, there a 2n-1 j states per n.
        matrixSize+=2*n-1;
        //States where mj>j are also non physical. This is not a problem if mj= 1/2 (notice that, in this case, the second
        //term is equal to 0). If, say, mj=3/2, we'd like to remove all j states where mj>j, or more specifically, where
        // j=1/2. There are 2 such states (^(1/2)S and ^(1/2)P). In this case, our subtracted term is equal to 2, as we
        // expect. This result applies generally.
        matrixSize-=(2*mj-1);
    }
    for(int n=nMax+1;n<=nMax+addedStates;n++)
    {
        //We wish to add the lower angular momentum states of higher n, as they may be relevant in our stark map. More
        //specifically, we add the S,P, and D states for higher n. For each n, there are 5 j states associated with these
        //angular momentum states. Again, we wish to subtract states where mj>j.
        matrixSize+= 5-(2*mj-1);
    }
    return matrixSize;
}

double StarkHamiltonian::ClebschGordanCoeff(int j1, int j2, int m1, int m2, int j, int m) {
	double coeff = 0;
	double factor = 0;
	int zmin, zmax;
	int par;
	double sum = 0;
	//printf("first: j2 = %d m2 = %d ((j2+m2)/2) = %d\n", j2, m2, ((j2 + m2) / 2));
	//Make sure the arguments are valid. If not, return zero.

	if ((((int) fabs(j1)) % 2 != ((int) fabs(m1)) % 2) || (((int) fabs(j2)) % 2 != ((int) fabs(m2)) % 2) || (((int) fabs(j)) % 2 != ((int) fabs(m)) % 2) || j < 0 || j1 < 0 || j2 < 0 || ((int) fabs(m)) > j || ((int) fabs(m1)) > j1 || ((int) fabs(m2)) > j2 || (j1 + j2) < j || ((int) fabs(j1 - j2))
			> j || (m1 + m2) != m)
		return coeff;
	//return 0;

	factor = binom(j1, (j1 + j2 - j) / 2) / binom((j1 + j2 + j + 2) / 2, (j1 + j2 - j) / 2);
	factor = factor * binom(j2, (j1 + j2 - j) / 2) / binom(j1, (j1 - m1) / 2);
	factor = (factor / binom(j2, (j2 - m2) / 2)) / binom(j, (j - m) / 2);
	factor = sqrt(factor);
	//printf("hi %f", factor);

	zmin = MAX(0, j2 + (j1 - m1) / 2 - (j1 + j2 + j) / 2);
	zmin = MAX(zmin, j1 + (j2 + m2) / 2 - (j1 + j2 + j) / 2);
	zmax = MIN((j1 + j2 - j) / 2, (j1 - m1) / 2);
	zmax = MIN(zmax, (j2 + m2) / 2);
	//printf("j2 = %d m2 = %d ((j2+m2)/2) = %d\n", j2, m2, ((j2 + m2) / 2));
	//printf("zmin = %d zmax = %d \n", zmin, zmax);
	//return 0;
	//<= or < ????
	for (int z = zmin; z <= zmax; z++) {
		//System.out.println("hi");
		par = 1;
		if (z % 2 != 0)
			par = -1;
		sum += par * binom((j1 + j2 - j) / 2, z) * binom((j1 - j2 + j) / 2, (j1 - m1) / 2 - z) * binom((-j1 + j2 + j) / 2, (j2 + m2) / 2 - z);
	}

	//printf(" sum=  %f", sum);

	coeff = factor * sum;
	//printf(" coeff=  %f\n", coeff);
	return coeff;
}

double StarkHamiltonian::binom(int n, int r) {
	if (n == r || r == 0)
		return 1;
	if (r == 1)
		return (double) n;
	return ((double) n / (double) (n - r)) * binom(n - 1, r);
}

double StarkHamiltonian::Wigner3jSymbol(int j1, int j2, int j3, int m1, int m2, int m3) {
	return (pow(-1, .5 * (j1 - j2 - m3)) * ClebschGordanCoeff(j1, j2, m1, m2, j3, -1 * m3)) / sqrt((double) (j3 + 1));
}

/**
 * Returns maximum of two integers
 */
int StarkHamiltonian::MAX(int n1, int n2) {
	if (n1 > n2)
		return n1;
	else
		return n2;
}

/**
 * Returns minimum of two integers
 */
int StarkHamiltonian::MIN(int n1, int n2) {
	if (n1 < n2)
		return n1;
	else
		return n2;
}

/**
 * returns 1 if equal, zero if not
 */
int StarkHamiltonian::delta(double a, double b){
        if(a==b)
            return 1;
        else
            return 0;
    }
