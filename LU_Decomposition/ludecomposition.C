//---------------------------
// We will write a program that performs LU decomposition on
// a predefined 5x5 matrix. We will write the original matrix, A.
// We will then write the LU decom matrix, L, U, LxU, the
// inverse of A=LU, and det(A).
//----------------------------
using namespace std;
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>

const int N=4;

void InitMat(double A[N][N], double B[N] ){
A[0][0] =   3.0;  A[0][1] =  4.0;  A[0][2] =  1.0;   A[0][3] = 3.0;  B[0] =   1;
A[1][0] =  -5.0;  A[1][1] =  4.0;  A[1][2] = -1.0;   A[1][3] = 4.0;  B[1] =   1;
A[2][0] =   8.0;  A[2][1] = -2.0;  A[2][2] =  6.0;   A[2][3] = 5.0;  B[2] =   1;
A[3][0] =  -4.0;  A[3][1] =  3.0;  A[3][2] = -0.5;   A[3][3] = 7.0;  B[3] =   1;
}

void PrintMatVec(double A[N][N], double B[N]){
	int i,j;

	for(i=0;i<N;i++){cout << setw(8) << "A[" << i << ",.] =" << setw(8) << " ";
        	for(j=0;j<N;j++){cout << setw(12) << A[i][j] << setw(12) << " ";}

		 cout << "    " << setw(10) << "B[" << i << "] =" <<  B[i] <<  endl;
	}
	}

void PrintMat(double A[N][N]){
        int i,j;

        for(i=0;i<N;i++){cout << setw(8) << "A[" << i << ",.] =" << setw(8) << " ";
                for(j=0;j<N;j++){cout << setw(12) << A[i][j] << setw(12) << " ";}
		cout << " " << endl;
        }
        }

void PrintMatInv(double I[N][N]){
        int i,j;

        for(i=0;i<N;i++){cout << setw(8) << "I[" << i << ",.] =" << setw(8) << " ";
                for(j=0;j<N;j++){cout << setw(12) << I[i][j] << setw(12) << " ";}
                cout << " " << endl;
        }
        }


void LUdec(double A[N][N]){
	int i,j,k;

	for(j=0;j<N;j++){
		for(i=0;i<j+1;i++){if(i>0){for(k=0;k<i;k++){A[i][j]-=A[i][k]*A[k][j];}}}
		if(j<N-1){for(i=j+1;i<N;i++){if(j>0){for(k=0;k<j;k++){A[i][j]-=A[i][k]*A[k][j];}}
			A[i][j]=A[i][j]/A[j][j];}}
	}
	}

void FSub(double L[N][N], double B[N], double x[N]){
	int i,j,k;

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
		if(i==j){L[i][j] = 1;}
		if(i<j){L[i][j] = 0;}
		}
	}
	cout << "For the Lower Triangular matrix L" << " " << endl;
	
        x[0]=B[0];
        for(k=0;k<N;k++){x[k]=B[k]; for(j=0;j<k;j++){x[k]=x[k]-L[k][j]*x[j];}}
	cout << " " << endl;
        for(i=0;i<N;i++){cout << setw(8) << "y[" << i << "] =" << " " << setw(8) << x[i] << endl;}
	}

void BSub(double U[N][N], double B[N], double x[N]){
	int i,j,k;

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){if(i>j){U[i][j] = 0;}}
	}
	cout << "For the upper triangular matrix U" << " " << endl;
	if(U[N-1][N-1]==0){cout << "error" << endl; }
	x[N-1]=B[N-1]/U[N-1][N-1];
	for(k=N-2;-k<1;k--){x[k]=B[k];for(j=k+1;j<N;j++){x[k]=x[k]-U[k][j]*x[j];}
      	  x[k]=x[k]/U[k][k];
	}
	cout << " " << endl;
	for(i=0;i<N;i++){cout << setw(8) << "x[" << i << "] =" << " " << setw(8) << x[i] << endl;}
	}

void test(double A[N][N], double B[N], double x[N], double o[N], double e[N]){
	int i,j,k;
	o[0]=0;

	InitMat(A,B);
	for(k=0;k<N;k++){for(i=0;i<N;i++){o[k]=o[k]+A[k][i]*x[i];}}
	for(k=0;k<N;k++){e[k]=B[k]-o[k];}

	cout << setw(8) << "A*x" << setw(12) << " " << "B" << " " << setw(15) << "error" << endl;
	cout << setw(8) << "---" << setw(12) << " " << "---" << " " << setw(15) << "------" << endl;

	for(j=0;j<N;j++){cout << setw(8) << o[j] << setw(12) <<  " " <<  B[j] << " " << setw(15) << e[j]  << endl;}
	
	cout << " " << endl;
	}

void LtimesU(double L[N][N], double U[N][N], double LU[N][N]){
	int i,j,k;

	for(i=0;i<N;i++){for(j=0;j<N;j++){for(k=0;k<N;k++){LU[i][j]+=L[i][k]*U[k][j];}}}
	cout << "L times U" << " " << endl;
	}

void inverse(double A[N][N], double L[N][N], double U[N][N], double B[N], double y[N], double x[N], double I[N][N], double ID[N][N]){
	int i,j,q,k;
	cout << " " << endl;
	cout << " " << endl;
	cout << "Calculating the Inverse of A" << endl;
	for(i=0;i<N;i++){
		InitMat(L,B);
		for(j=0;j<N;j++){if(i==j){B[j]=1;}else{B[j]=0;}}
		cout << " " << endl;
		cout << setw(12) << "For B = " << " ";
		for(q=0;q<N;q++){cout << setw(6)  << B[q] << " ";}
		cout << " " << endl;
		LUdec(L);
		FSub(L,B,y);
		InitMat(U,B);
		LUdec(U);
		BSub(U,y,x);
		for(j=0;j<N;j++){I[j][i]=x[j];}
	}
	cout << "Inverse of A" << " " << endl;
	PrintMatInv(I);
	cout<< " " << endl;
	cout << "Test to see if AA^(-1) = 1." << " " << endl;
	InitMat(A,B);
	for(k=0;k<N;k++){
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){ID[i][k]+= A[i][j]*I[j][k];
	}
	}
	}
	for(i=0;i<N;i++){for(j=0;j<N;j++){if(ID[i][j]<1e-10){ID[i][j]=0;}}}
	cout << "A times its Inverse Yields" << " " << endl;
	PrintMat(ID);
	}

void det(double U[N][N], double B[N], double y[N], double x[N]){
	int i;
	double p;
	p=1;

	cout << " " << endl;
        cout << " " << endl;
	cout << "Calculating the determinant of A" << " " << endl;
	cout << " " << endl;
	InitMat(U,B);
	LUdec(U);
	BSub(U,y,x);
	for(i=0;i<N;i++){p*=U[i][i];}
	cout << " " << endl;
	cout << "The determinant of the matrix A is " << p << ".  " << endl;
	}



int main(){

double A[N][N], U[N][N], L[N][N], LU[N][N], B[N], y[N], x[N], o[N], e[N], I[N][N], ID[N][N];
int i,j,k;
double M;

InitMat(L,B);
PrintMat(L);
cout << "LU Decomposition" << " " << endl;
LUdec(L);
PrintMat(L);
FSub(L,B,y);
cout << " " << endl;
PrintMatVec(L,B);
cout << " " << endl;
InitMat(U,B);
LUdec(U);
BSub(U,y,x);
cout << " " << endl;
//test(A,B,x,o,e);
PrintMatVec(U,B);
cout << " " << endl;
LtimesU(L,U,LU);
PrintMat(LU);
inverse(A,L,U,B,y,x,I,ID);
det(U,B,y,x);
}


