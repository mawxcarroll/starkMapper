#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include "mpi.h"
#include "netcdf.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "RydbergAtom.h"
#include "StarkHamiltonian.h"

#define DEBUG 0
#define DATAOUT 1
#define STATESOUT 0

#define NTYPES 3
#define PUATOM 0
#define PLATOM 1
#define SATOM  2

#define ALLPU 1
#define SUPERPOSITION 2

using namespace std;


//handles pnetcdf file errors
static void handle_error(int status);



extern "C" {
double pdlamch_(int* context, char* c);
double pdelget_(char* scope, char* top, double* alpha, double* A,int* i, int* j, int* descA);
void pdelset_(double* d1, int* i1, int* i2, int* i3, double* d2);
void pdsyevx_(char *, char *, char *, int *, double *, int *, int *, int *, double *, double *, int *, int *, double *, int *, int *, double *, double *, double *, int *, int *, int *, double *, int *, int *, int *, int *, int *, double *, int *);
void pdtran_(int * M, int * 	N, double * ALPHA, double * A, int * IA, int * JA, int * DESCA, double * BETA, double * C, int * IC, int * JC, int * DESCC);
void pdcopy_(int * N, double * X, int * IX, int * JX, int * DESCX, int * INCX, double * Y, int * IY, int * JY, int * DESCY, int * INCY);
void pdscal_(int * N, double * ALPHA, double * X, int * IX, int * 	JX, int * DESCX, int * INCX);
void pdgemm_(char * TRANSA, char * TRANSB, int * M,int * N, int * K, double * ALPHA, double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB, double * BETA, double * C, int * IC, int * JC, int * DESCC);
void pdgeadd_(char *	TRANS, int * M, int * N, double * ALPHA, double * A, int * IA, int * JA, int * DESCA, double * BETA, double * C, int * IC, int * JC, int * DESCC, int sLen);
void pdlaprnt_(int *, int *, double *, int *, int *, int *, int *, int *, char *, int *, double *, int slen);

int numroc_(int* i1, int* i2, int* i3, int* i4, int* i5);
void descinit_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
void descset_(int *, int *, int *, int *, int *, int *, int *, int *, int *);

void Cblacs_gridexit (int ConTxt);
void Cblacs_gridinfo(int , int *, int *, int *, int *);
void Cblacs_gridinit(int *ConTxt, char *order, int nprow, int npcol);
void Cblacs_pinfo(int *mypnum, int *nprocs);
void Cblacs_get(int context, int what, int *val);
void Cblacs_gridinit(int* context, char* order,int nproc_rows, int nproc_cols);
void Cblacs_pcoord(int context, int p, int* my_proc_row, int* my_proc_col);
void Cblacs_exit(int doneflag);
}
//extern void pslaprnt_(int *m, int *n, double *a, int *ia, int *ja, int *descrip, int *code1, int *code2, char *mname, size_t *mlength, int *code3, double *wrk);

int main(int argc, char **argv) {
	cout << "In c++ file" << endl;
	int i, j, k;
	double tconv = 4.134137174427E10;
	int retval, ncid, varid;


	/************  MPI ***************************/
	int myrank_mpi, nprocs_mpi;
	MPI_Init( &argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
	/*********************************************/

	int c = 1;
	/*** Command Line Arguments ******************/
	//name of the save file
	char *saveFile = argv[c++];
	//atom file
	char *atomFile = argv[c++];
	//minimum value of n
	int nMin = atoi(argv[c++]);
	//maximum value of n
	int nMax = atoi(argv[c++]);
	//mj
	double mj = atof(argv[c++]);
	//electric field in V/cm
	double efieldMin = atof(argv[c++]);
	double efieldMax = atof(argv[c++]);
	double efieldRes = atof(argv[c++]);
	string saveDir = string(argv[c++]);
    if(saveDir[saveDir.size()-1] != '/')
        saveDir.append("/");


	if(DEBUG==1) {
		printf("nMin = %d\n", nMin);
		printf("nMax = %d\n", nMax);
		printf("mj = %f\n", mj);
		printf("savefile = %s\n", saveFile);
		printf("atomFile = %s\n", atomFile);
	}

	//Rectangle rect(5,4);
	//cout << "******************  " << rect.area() << "\n";

	//timing stuff
	struct timeval t1;
	gettimeofday(&t1, NULL);

	RydbergAtom ra(atomFile); //create the Rydberg atom object, which stores polarizability and quantum defects

    //cout << "Polarizability = " << ra.polarizability << "\n";
    //create the class that calculates matrix elements for the Hamiltonian
	StarkHamiltonian sh(0.01, pow(9.023,(1./3.)), 10E-10, 10E-5, nMin, nMax, mj,3);

    int N = sh.getSizeOfMatrix();

    //cout << N << "\n";
    //cout << sh.getMatrixElement(1, 15, 0, 14, 0.5, 1.5, 0.5, 2.6417,3.131 , 0.0) << "\n";
    //cout << sh.getMatrixElement(2, 16, 1, 16, 0.5, 1.5, 1.5, 1.346, 2.6417, 1.0) << "\n";
    //cout << sh.getMatrixElement(2, 16, 1, 16, 0.5, 2.5, 1.5, 1.346, 2.6417, 1.0/sh.eFieldConversion) << "\n";
    //cout << sh.getMatrixElement(12, 16, 11, 15, 0.5, 12.5, 11.5, 1.346, 2.6417, 1.0/sh.eFieldConversion) << "\n";



	/************  BLACS *************************/
	int ictxt, nprow, npcol, myrow, mycol,nb;
	int info,itemp;
	int ZERO=0,ONE=1,SIX=6;
	double MONE = -1.0, FONE=1.0, FZERO=0.0;

	//define grid
	//size of grid depends on size of matrix
	//assume small matrix size for now

	nprow = 1;
	npcol = 1;


	//printf("%d, %d\n",nprow,npcol);
	nb = 128;


	Cblacs_pinfo( &myrank_mpi, &nprocs_mpi ) ;
	Cblacs_get( -1, 0, &ictxt );
	Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

	/*********************************************/




    /************* Allocate memory for diagonalization and time evolution ***********************/

	int descA[9],descZ[9], descP[9], descQ[9], descIR[9], descII[9], descB[9];
	int mA = numroc_( &N, &nb, &myrow, &ZERO, &nprow );
	int nA = numroc_( &N, &nb, &mycol, &ZERO, &npcol );
	int nx = numroc_( &N, &nb, &myrow, &ZERO, &nprow );
	int my = numroc_( &N, &nb, &myrow, &ZERO, &nprow );
	//printf("mA = %d and nA = %d\n", mA,nA);
	//matrix A is the Hamiltonian
	descinit_(descA, &N,   &N,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA,  &info);
	//matrix Z holds the eigenvectors
	descinit_(descZ, &N,   &N,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA,  &info);
	//matrix Q for calculation of U(t)
	/**
	 * How much of Q we need depends on the initial state.  If the initial state
	 * has n nonzero entries, then we need n columns of Q. For this calculation,
	 * we'll generally need all of Q. The state vector will evolve to have multiple
	 * nonzero entries and we are not free to reorder them.
	 */
    int nInitial = N;
	descinit_(descQ, &N,   &nInitial,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA,  &info);
	descinit_(descB, &N,   &nInitial,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA,  &info);
	//column vectors P and I
	descinit_(descP, &N,   &ONE,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA,  &info);
	descinit_(descIR, &nInitial,   &ONE,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA,  &info);
	descinit_(descII, &nInitial,   &ONE,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA,  &info);

	double *A = (double*) malloc(mA*nA*sizeof(double));
	double *B = (double*) malloc(mA*nA*sizeof(double));
	double *Z = (double*) malloc(mA*nA*sizeof(double));
	double *P = (double*) malloc(mA*sizeof(double));
	double *II = (double*) malloc(mA*sizeof(double));
	double *IR = (double*) malloc(mA*sizeof(double));
	double *Q = (double*) malloc(mA*nA*sizeof(double));
	double *W = (double*) malloc(N*sizeof(double));
	double *psiR = (double*) malloc(N*sizeof(double));
	double *psiI = (double*) malloc(N*sizeof(double));

	//define the sizes of the fortran workspaces
	int LWORK=1;
	int LIWORK=1;
	int *IWORK = (int*) malloc(LIWORK*sizeof(int));
	double *WORK = (double*)  malloc(LWORK*sizeof(double));
	int *IFAIL = (int*)  malloc(N*sizeof(int));
	int *ICLUSTR = (int*) malloc(N*sizeof(int));
	double *GAP = (double*) malloc(N*sizeof(double));

	ofstream myFile;
    myFile.open ("MTDTest.mtd");
    if (myFile.is_open())
    {
        myFile << N;
        myFile.close();
    }


	/*************************************************************************************/

	/************************** Create the Hamiltonian ***********************************/

	/**
	 * We're using ScaLAPACK to do the numerical heavy-lifting. In the memory allocation
	 * section, we've already set up all of the memory we need to do a full scale time
	 * evolution problem. To create a Stark map, we just need to call pdsyevx to find
	 * the eigenvalues (it also finds the eigenvectors, but we don't really need that for
	 * Stark maps).
	 *
	 * The memory we'll be using to store the Hamiltonian matrix is addressed by the
	 * pointer A. We need to loop over the elements of the Hamiltonian in some way similar
	 * to the java code from the previous version. Many matrix elements we can just set
	 * to zero, the rest we can calculate by calling the getMatrixElement() method of
	 * the StarkHamiltonian class.
	 *
	 * We will set the matrix elements using the ScaLAPACK function pdelset(),
	 *      pdelset_(A, &i, &j, descA, &element)
	 * where
	 *      A = our Hamiltonian matrix
	 *      i = row index
	 *      j = column index
	 *      descA = matrix descriptor (contains information about how the memory is allocated)
	 *      element = the matrix element.
	 * Note that this is calling, at root, a Fortran function. Thus, we need to pass everything
	 * by reference (hence the &; A and descA are already pointers). Also, Fortran indexes
	 * arrays from 1, not 0 like C/C++/Java. So the upper left-most element of the matrix
	 * is the (1,1) element NOT the (0,0) element. Finally, since the matrix is symmetric, we
	 * can immediately call
	 *      pdelset_(A, &j, &i, descA, &element)
	 * and our loop only needs to loop over the upper right triangle of the matrix (and the
	 * diagonal).
	 */
    for(double ee = efieldMin; ee < efieldMax; ee+= efieldRes)
    {
        // nStart is the row of the matrix that the current n value starts at - 1.
        int nStart=0;
        // F is the field
        double F=ee/sh.eFieldConversion;
        // Loop over range of n, this is used to recover the principle quantum number that we are interested in
        for(int n=nMin;n<=nMax+sh.addedStates;n++)
        {
           // Loop over all possible j values for a given n
           int jMax;
           //if its not an added state (equivalently, if all j are included)
           if (n<=nMax)
               jMax = 2*n-1-(int)(2*mj-1);
           //if it is, we only want to account for the lower l states (all physical s,p,d)
           else
           // There are a maximum of 5 j-states for minimal mj (mj = 0.5) for the s, p, d orbitals;
           // we subtract off the non-physical j-states for a given mj
               jMax= 5- (2*mj-1);
           for(int jI=1;jI<=jMax;jI++)
           {
               //::: FORMULAS TO FIND L's and J's :::
               //        m_j=3/2
               // jI        L          J
               // 1         1          1.5
               // 2         2          1.5
               // 3         2          2.5
               // 4         3          2.5
               // 5         4          3.5
               // ..        ..         ..
               //We observe that l=floor(jI/2)+floor(mj). Notice that integer division and casting to an int perform the floor functions for us.
               //Additionally, j = jI/2 -(1+(-1)^jI)/4+floor(mj). For odd jI, the second term is 0, for even jI, it is 1/2.
               //These formulas hold generally for all mj.
               // Use the jIndex and mj to recover the orbital angular momentum quantum number
               int l= jI/2 + (int)mj;
               // Use the jIndex to recover the total angular momentum quantum number
               double j= jI/2.0 - (1+pow(-1,jI))/4.0 +(int) mj;
               // find the row index... there is more than one state for a given j, so jIndex is necessary
               int rowIndex=nStart+jI;
               //analagous to nStart, but used for np instead.
               int npStart=0;
               // We utilize the Hermiticity of the Hamiltonian and simply find the components of H in the lower left triangle.
               for(int np=nMin;np<=n;np++)
               {
                   //if the n's are the same, then we do not need to go over all j states; we need only loop up until the jIth column in that n.
                   //ADD IN JPMAX TO TAKE CARE OF ADDED STATES
                   int jpMax;
                   if(np<=nMax)
                       jpMax=2*np-1-(int)(2*mj-1);
                   else
                       jpMax=5-(2*mj-1);
                   if(np==n)
                   {
                       for(int jpI=1;jpI<=jI;jpI++)
                       {
                           //analogous to formulas for l and j
                           int lp= jpI/2 + (int) mj;
                           double jp= jpI/2.0 - (1+pow(-1,jpI))/4.0 + (int) mj;
                           int colIndex=npStart+jpI;
                           //calculates the matrix element
                           double element = sh.getMatrixElement(double(l),double(n),double(lp),double(np),mj,j,jp,ra.getQuantumDefect(n,l,j),ra.getQuantumDefect(np,lp,jp),F);
                           //if ((l!=lp+1 || l!=lp-1) && rowIndex!=colIndex)
                           //    element=0;
                           //if(element!=0)
                             //  cout << "{"<<rowIndex<<","<<colIndex<<","<< n << "," << l <<","<<j<<","<<ra.getQuantumDefect(l,j)<<","<<np<<","<<lp<<","<<jp<<","<<ra.getQuantumDefect(lp,jp)<<","<<element<<"}"<<"\n";
                           // if(rowIndex!=colIndex)
                            //   element=0;
                           //cout << n << ", " << l << ", " << j << ", " << np << ", " << lp << ", " << jp << "--> " << element << "\n";
                           // The inner product is real for all states, so the matrix is symmetric
                           //sets the matrix element to the (i,j) entry and the (j,i) entry
                           pdelset_(A,&rowIndex,&colIndex,descA,&element);
                           pdelset_(A,&colIndex,&rowIndex,descA,&element);
                       }
                   }
                   //if they aren't the same, then we need to do the same as before, but over all j values in np.
                   else
                   {
                       for(int jpI=1;jpI<=jpMax;jpI++)
                       {
                           int lp= jpI/2 + (int) mj;
                           double jp= jpI/2.0 - (1+pow(-1,jpI))/4.0 + (int) mj;
                           int colIndex=npStart+jpI;
                           double element = sh.getMatrixElement(double(l),double(n),double(lp),double(np),mj,j,jp,ra.getQuantumDefect(n,l,j),ra.getQuantumDefect(np,lp,jp),F);
                           //if(element!=0)
                               //cout << "{"<<rowIndex<<","<<colIndex<<"," << n << "," << l <<","<<j<<","<<ra.getQuantumDefect(l,j)<<","<<np<<","<<lp<<","<<jp<<","<<ra.getQuantumDefect(lp,jp)<<","<<element<<"}"<<"\n";
                           //cout << "["<< rowIndex << ", " <<  colIndex << "] --> "<< n << ", " << l << ", " << j << "," << ra.getQuantumDefect(l,j) << "  |||  " << np << ", " << lp << ", " << jp << ", " << ra.getQuantumDefect(lp,jp) <<" ||||| "<<element<< "\n";
                           //if(rowIndex!=colIndex)
                            //  element=0;
                           //cout << n << ", " << l << ", " << j << ", " << np << ", " << lp << ", " << jp << "--> " << element << "\n";
                           // The inner product is real for all states, so the matrix is symmetric
                           pdelset_(A,&rowIndex,&colIndex,descA,&element);
                           pdelset_(A,&colIndex,&rowIndex,descA,&element);
                           // This is for writing the matrix to a file... Fortran starts @ 1, while C++ indexes from 0
                       }
                   }
                   npStart+=jpMax;
               }
           }
           nStart+=jMax;
        }

        /*************************************************************************************/
        //timing stuff
        struct timeval t2;
        gettimeofday(&t2, NULL);
        //printf("time elapsed: %d secs\n",(t2.tv_sec-t1.tv_sec));

        /**
         *  Find the eigenvalues. The solver, pdsyevx, is called twice. The first call
         *  returns optimal sizes for the workspaces. The second call solves the system.
         *  For generating a Stark Map we should be able to use a solver that does not
         *  save the eigenvectors. For the current LZ problem, we'll need the eigenvectors.
         */


        char *matrixname="test";
        char cmatnm = 'A';
        int len_m = strlen(matrixname);
        //pdlaprnt_(&N, &N, A, &ONE, &ONE, descA, &ZERO, &ZERO, matrixname, &SIX, WORK,len_m);
        double alpha = 1.0; double beta = 0.0;
        int num13 = 13;
        int negnum13=-13;
        double ABSTOL=pdlamch_(&ictxt, "U");
        //printf("ABSTOL = %.32f\n", ABSTOL);
        int MM, NZ; //what are these for?
        int INFO;
        //printf("HI!\n");
        /**********************************
         * set LIWORK and LWORK to -1     *
         * when pssyevx is called, it will*
         * return the optimal values in   *
         * WORK[0] and IWORK[0]
         **********************************/

        LIWORK=-1;
        LWORK=-1;
        char vv = 'V';
        char aa = 'A';
        char uchar = 'U';
        char ss = 'S';

        pdsyevx_(&vv,&aa,&uchar,&N,A,&ONE,&ONE,descA,&FZERO,&FZERO,&ZERO,&ZERO,&ABSTOL, &MM, &NZ,W,&MONE,Z,&ONE,&ONE,descZ,WORK,&LWORK,IWORK,&LIWORK, IFAIL, ICLUSTR, GAP, &INFO);
        LIWORK = IWORK[0]+1;
        LWORK = (int)(WORK[0])+1;//+7964*N;
        //printf("LIWORK=%d \n",LIWORK);
        //printf("LWORK=%d \n",LWORK);
        //LWORK *= 4;
        //LIWORK *= 4;
        IWORK = (int*) malloc(LIWORK*sizeof(int));
        WORK = (double*)  malloc(LWORK*sizeof(double));
        IFAIL = (int*)  malloc(N*sizeof(int));
        ICLUSTR = (int*) malloc(N*sizeof(int));
        GAP = (double*) malloc(N*sizeof(double));

        pdsyevx_(&vv,&aa,&uchar,&N,A,&ONE,&ONE,descA,&FZERO,&FZERO,&ZERO,&ZERO,&ABSTOL, &MM, &NZ,W,&MONE,Z,&ONE,&ONE,descZ,WORK,&LWORK,IWORK,&LIWORK, IFAIL, ICLUSTR, GAP, &INFO);
        //for(int i=0;i<N;i++)
        //    oneFieldWithEigenvalues.push_back(W[i]*sh.energyConversion);

        /********************************************** Write To File ***********************************************/

        stringstream sss, nMinStrS, nMaxStrS, mjStrS;
        sss << ee;
        string fieldString = sss.str();
	nMinStrS << nMin;
	string nMinString = nMinStrS.str();
	nMaxStrS << nMax;
	string nMaxString = nMaxStrS.str();
	mjStrS << mj;
	string mjString = mjStrS.str();

        ofstream outFile;
        string smSaveFile = saveDir+"Stark_Map_n_"+nMinString+"_"+nMaxString+"_mj_"+mjString+"_F_"+fieldString+".dat";
        outFile.open(smSaveFile.c_str(), ios::out | ios::binary);
        outFile.write((char*)&F, sizeof(double)); //prints field value as the first entry in the file
        for(i=0;i<N;i++)
        {
            double convEigenvalues=W[i]*sh.energyConversion;
            outFile.write((char*)&convEigenvalues,sizeof(double)); //prints out Eigenvalues horizontally
        }
	outFile.close();

	//Write out Full Hamiltonian
	/*
        string hSF = saveDir+"Hamiltonian_"+nMinString+"_"+nMaxString+"_mj_"+mjString+"_F_"+fieldString+".dat";
        outFile.open(hSF.c_str(), ios::out | ios::binary);
        outFile.write((char*)&F, sizeof(double)); //prints field value as the first entry in the file
        char defaultTop=' ';
        char updateAllProcesses='A';
        for(i=1; i<=N; i++)
	    for(j=1;j<=N;j++)
	    {
		    double el=pdelget_(&updateAllProcesses,&defaultTop,&el,A,&i,&j,descA);
		    outFile.write((char*)&el,sizeof(double)); //prints element to file
	    }
        outFile.close();
        */
	cout << "The file " << smSaveFile << " Has been successfully written" << endl;
    }
	Cblacs_gridexit( 0 );
	MPI_Finalize();
	free(A);
	free(W);
	free(Z);
	free(IWORK);
	free(WORK);
	free(IFAIL);
	free(ICLUSTR);
	free(GAP);
	return 0;
}

static void handle_error(int status)
{
	printf("%s\n", nc_strerror(status));
	//exit(0);
}

