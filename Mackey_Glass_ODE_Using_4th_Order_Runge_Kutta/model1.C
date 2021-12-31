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

const int resSize=400;
const int inSize=1;
const int tfin=1000;
const double alp=0.3;
const double wspec=1.3;
double W[resSize][resSize];
double A[resSize][resSize];

//ofstream prt1("model1.dat");
ofstream prt1("model1scale13.dat");
ofstream prteigs("model1eigs.dat");

int main(){
	double xnew[resSize],xold[resSize],maxlam,s;
	int i,j,k,t,info;
	int nn=resSize;
	char jobvl='N', jobvr='N';
	const int lwork=8*resSize;
	double valr[resSize],vali[resSize],work[lwork],dumv[2];

	t=0;
	srand48(1234);
	for(i=0;i<resSize;i++){
		xold[i]=(drand48()-0.5)*2.0; 
		for(j=0;j<resSize;j++){
			W[i][j] = (drand48()-0.5)*2.0; 
		}
	}

	for(i=0;i<resSize;i++){ for(j=0;j<resSize;j++){ A[i][j]=W[i][j];}}

	info=LAPACKE_dgeev_work(LAPACK_ROW_MAJOR,jobvl,jobvr,nn,A[0],nn,valr,vali,dumv,nn,dumv,nn,work,lwork);
	maxlam=0.0;
	for(i=0;i<resSize;i++){
		s=sqrt(valr[i]*valr[i]+vali[i]*vali[i]);
		if(s>maxlam){ maxlam=s; }
		prteigs << valr[i] << " " << vali[i] << endl; 
	}

	for(i=0;i<resSize;i++){ for(j=0;j<resSize;j++){ W[i][j]=wspec*W[i][j]/maxlam; }}

	while(t<tfin){

	for(i=0;i<6;i++){ prt1  << xold[i] << " "; } prt1 << endl;

	for(i=0;i<resSize;i++){
		xnew[i]=0.00;
		for(j=0;j<resSize;j++){
			xnew[i]+=W[i][j]*xold[j];
		}
	        xnew[i]=(1.0-alp)*xold[i]+alp*tanh(xnew[i]);
	}
	
	t++;
	for(i=0;i<resSize;i++){ xold[i]=xnew[i]; } 
	} //end of while loop



}
