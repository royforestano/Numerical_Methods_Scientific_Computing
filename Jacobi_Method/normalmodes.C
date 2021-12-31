/*--------------------------------------------------------------
 ODE: Solution
-------------------------------------------------------------- */

using namespace std;
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <mkl_lapacke.h>
#include <fstream>
#include <ctime>

const int N=8;
clock_t begintime,endtime;
double elas_time;
//double rowmodify = N/2;
//double r = rowmodify;
double PI = 3141592653589793238;

void InitMat(double A[N][N]){
	int i,j,left,right;
	for(i=0;i<N;i++){
		if(i>0){left=i-1;}else{left=N-1;}
		if(i<N-1){right=i+1;}else{right=0;}
		for(j=0;j<N;j++){
			if(i==j){
				A[i][j]=2.0;
			}else if(j==left){
				A[i][j]=-1.0;
			}else if(j==right){
				A[i][j]=-1.0;
			}else{
				A[i][j]=0.0;
	}}}

	A[N/2][N/2]=A[N/2][N/2]/0.01;
}


void PrintMat(double A[N][N]){
	int i,j,k;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
		cout << setw(14) << A[i][j] << " ";
		}
	cout << endl;
	}}

int main(){
	double A[N][N],val[N];
	int i,j,k,t,info;
	int nn=N;
	const int lwork=8*N;
	char jobvl='N', jobvr='N';
	double  work[lwork];

	InitMat(A);
	cout << " " << "This is the linear equation relating amplitudes of normal modes" << " " << endl;
	PrintMat(A);

// For determining the time it takes to diagonalize the matrix ------------------------------------------------
	cout << endl;
	begintime = clock();
	info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',nn,A[0],nn,val);
	endtime=clock();


//restore units to frequency and set a to unit length-------------------------------------------------- 
	for(i=0;i<N;i++){
		val[i]=2*(1-cos(2*PI*i/N))*val[i];
	}

// Printing the eigenvalues -----------------------------------------------------------------------------------
	cout << " " << "Here are the eigenvalues" << " " << endl;
	for(i=0;i<N;i++){
		cout << " " << setw(8) << "w" << i+1 << "^2 =" << " " << setw(12) << val[i] << " " << "s^-2" << endl;
	}

// Printing the eigenvectors -----------------------------------------------------------------------------------
	cout << endl;
	cout << "Each column is the eigenvector corresponding to the eigenvalues above, respctively" << " " << endl;
	PrintMat(A);

// Printing the time it takes to diagonlaize the matrix --------------------------------------------------------	
	cout << endl;
	cout << "Time Spent" << " " << endl;
	elas_time = (double)(endtime-begintime)/CLOCKS_PER_SEC;
	cout << "time in seconds: " << elas_time << endl;

}

