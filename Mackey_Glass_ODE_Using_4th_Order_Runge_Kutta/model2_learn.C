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

const int resSize=500;
const int inSize=1;
const int tfin=1000;
const double alp=0.3;
const double wspec=1.3;
const double reg=1.0E-8;
double W[resSize][resSize];
double newW[resSize][resSize];
double A[resSize][resSize];
double Win[resSize][inSize];
const int outSize=inSize;
double Wout[outSize][resSize];
const int transientTimes=100; // or 200
const int sampleTimes=500; 
double XX[resSize][sampleTimes];
double XXT[sampleTimes][resSize];
double XXXXT[resSize][resSize];
double INV[resSize][resSize];
double YY[outSize][sampleTimes];
double YYXXT[outSize][resSize];
double WoutX[outSize];
double WA[resSize][resSize];

double abs(double x){ return sqrt(x*x); }

ifstream readin("mglass.dat");
ofstream prtx("model2_learn.dat");
ofstream prtvalid("model2_learn_valid.dat");
ofstream prteigs("model2_learn_eigs.dat");


void MultMatVec1(double W[resSize][resSize], double xold[resSize], double xnew[resSize]){
	int i,j;
	for(i=0;i<resSize;i++){
		xnew[i]=0.0;
		for(j=0;j<resSize;j++){xnew[i]+=W[i][j]*xold[j];	}}}

void MultMatVec2(double Win[resSize][inSize], double u[inSize], double x[resSize]){
	int i,j;
	for(i=0;i<resSize;i++){
		x[i]=0.0;
		for(j=0;j<inSize;j++){x[i]+=Win[i][j]*u[j];	}}}

void MultMatMat(double A[resSize][resSize], double B[resSize][resSize], double C[resSize][resSize]){
        int i,j,k;
        for(i=0;i<resSize;i++){
                for(j=0;j<inSize;j++){
			C[i][j]=0.0;
			for(k=0;k<inSize;k++){
				C[i][j]+=A[i][k]*B[k][j];	}}}}

void Init(double x[resSize]){
        int i,j; double a,b,c;
        srand48(1234);
        for(i=0;i<resSize;i++){
                x[i]=(drand48()-0.5)*0.4;
                for(j=0;j<resSize;j++){
                        W[i][j] = (drand48()-0.5)*2.0;
                }
                for(j=0;j<inSize;j++){
                        Win[i][j] = (drand48()-0.5)*2.0;
                }
        }
        //for(i=0;i<400;i++){readin >> a  >> b  >> c ;}
}


void rescaleW(){
	double maxlam,s;
	int i,j,info;
	int nn=resSize;
	const int lwork=8*resSize;
	char jobvl='N', jobvr='N';
	double valr[resSize], vali[resSize], work[lwork], dumv[2];

	for(i=0;i<resSize;i++){for(j<0;j<resSize;j++){ A[i][j]=W[i][j]; }}
	info = LAPACKE_dgeev_work(LAPACK_ROW_MAJOR,jobvl,jobvr,nn,A[0],nn,valr,vali,dumv,nn,dumv,nn,work,lwork);
	maxlam=0.0;
	for(i=0;i<resSize;i++){
		s=sqrt(valr[i]*valr[i]+vali[i]*vali[i]);
		if(s>maxlam){maxlam=s;}
			prteigs << valr[i] << " " << vali[i] << " " << endl;
	}
	for(i=0;i<resSize;i++){ for(j=0;j<resSize;j++){ W[i][j] = wspec*W[i][j]/maxlam; }}
	//cout << "maxlam = " << maxlam << " expected " << sqrt(double(resSize))*0.558 << endl;
}

void makeWA(double A[resSize][resSize],double B[resSize][inSize],double C[outSize][resSize]){
        int i,j,k;
        double WinWout[resSize][resSize];

// Win times Wout ----------------------------------------------------------------------------
        for(i=0;i<resSize;i++){
                for(j=0;j<resSize;j++){
                        for(k=0;k<inSize;k++){
                                WinWout[i][j]+=B[i][k]*C[k][j];         }}}


        for(i=0;i<resSize;i++){
                for(j=0;j<resSize;j++){
                        WA[i][j]=A[i][j]+WinWout[i][j];                 }}}


int main(){
	double xnew[resSize],xold[resSize],xprime[resSize],maxlam,s,u[inSize], tin, uin, uintau;
	double xcopy[resSize], diff;
	int i,j,k,p,t,info;
	const int lwork=8*resSize;
	char jobvl='N', jobvr='N';
	double valr[resSize],vali[resSize],work[lwork],dumv[2];
	int pivotArray[resSize];

	t=0;
	Init(xold);
	rescaleW();

// loop time transient -----------------------------------------------------------------
	while(t<transientTimes){
	for(i=0;i<6;i++){ prtx << xold[i] << " "; } prtx << endl;

	readin >> tin >> u[0] >> uintau;

        MultMatVec1(W,xold,xnew);
	MultMatVec2(Win,u,xprime);

	for(i=0;i<resSize;i++){
	        xnew[i]=(1.0-alp)*xold[i]+alp*tanh(xnew[i]+xprime[i]);
		xold[i]=xnew[i];
	}

	t++;	
	} 	// This ends the while loop for loop time transient
	
	for(i=0;i<resSize;i++){ xcopy[i]=xold[i]; }

// loop time sample -------------------------------------------------------------------------------
        while(t<transientTimes+sampleTimes){
        for(i=0;i<6;i++){ prtx << xold[i] << " "; } prtx << endl;

        readin >> tin >> u[0] >> uintau;

        MultMatVec1(W,xold,xnew);
        MultMatVec2(Win,u,xprime);

        for(i=0;i<resSize;i++){
                xnew[i]=(1.0-alp)*xold[i]+alp*tanh(xnew[i]+xprime[i]);
		xold[i]=xnew[i];
		XX[i][t-transientTimes]=xnew[i];
        }
        
	for(i=0;i<outSize;i++){ YY[i][t-transientTimes]=u[i]; }

	t++;
	} 	// This ends the while loop for loop time sample

// construct Wout -------------------------------------------------------------------------------------
	for(i=0;i<resSize;i++){
		for(p=0;p<sampleTimes;p++){
			XXT[p][i]=XX[i][p];	}}
	
	for(i=0;i<resSize;i++){
		for(j=0;j<resSize;j++){
			XXXXT[i][j]=0.0;
			for(p=0;p<sampleTimes;p++){
				XXXXT[i][j]+=XX[i][p]*XXT[p][j];
	}}}
	
	for(i=0;i<resSize;i++){ XXXXT[i][i]+=reg;}

// constructing the Inverse of XXXXT-reg -----------------------------------------------------------------
	for(i=0;i<resSize;i++){ for(j=0;j<resSize;j++){ INV[i][j]=XXXXT[i][j]; }}
	info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR,resSize,resSize,INV[0],resSize,pivotArray);
	info=LAPACKE_dgetri_work(LAPACK_ROW_MAJOR,resSize,INV[0],resSize,pivotArray,work,lwork);
	MultMatMat(INV,XXXXT,A);
	cout << " A[i][i]= "; for(j=0;j<5;j++){ for(i=0;i<5;i++){ cout << A[i][i] << " "; } cout << endl;}
	cout << "A[i][i+1]= "; for(i=0;i<6;i++){ cout << A[i][i+1] << " "; } cout << endl;

// constructing YXT---------------------------------------------------------------------------------------
	for(i=0;i<outSize;i++){
                for(j=0;j<resSize;j++){
                        YYXXT[i][j]=0.0;
                        for(p=0;p<sampleTimes;p++){
                                YYXXT[i][j]+= YY[i][p]*XXT[p][j];
        }}}

// multiplying matrix elements of YXT on the left and INV on right----------------------------------------
	for(i=0;i<outSize;i++){
		for(j=0;j<resSize;j++){
			Wout[i][j]=0.0;
			for(p=0;p<sampleTimes;p++){
				Wout[i][j]+= YYXXT[i][p]*INV[p][j];
	}}}	
	cout << "Wout[0][i]= "; for(i=0;i<6;i++){ cout << Wout[0][i] << " ";} cout << endl;
	
// constructing Wout*x -----------------------------------------------------------------------------------
	t=transientTimes;
	diff=0.0;
	while(t<transientTimes+900){
		for(i=0;i<resSize;i++){ xold[i]=XX[i][t-transientTimes]; }
		cout << "Test Wout " << t << " ";
		s = 0.0;
		for(i=0;i<outSize;i++){
			u[i]=0.0;
			for(j=0;j<resSize;j++){
				WoutX[i]+=Wout[i][j]*xold[j];
			}
		prtvalid << t  << " " << WoutX[i] << " ";
		prtvalid << YY[i][t-transientTimes] << " " << abs(YY[i][t-transientTimes]-WoutX[i]) << endl;
		s=(YY[i][t-transientTimes]-WoutX[i])*(YY[i][t-transientTimes]-WoutX[i]);	
		}
                diff=sqrt(s);
                cout   << " " << diff << endl;
		t++;

	}

	//cout << "XX[0][i]= "; for(i=0;i<6;i++){ cout << XX[0][i] << " "; } cout << endl;
        //cout << "XX[1][i]= "; for(i=0;i<6;i++){ cout << XXT[1][i] << " "; } cout << endl;
        //cout << "XX[resSize-1][i]= "; for(i=0;i<6;i++){ cout << XXT[resSize-1][i] << " "; } cout << endl;

}
