/*--------------------------------------------------------------
 ODE: Solution
        Model:  QUANTUM Simple Harmonic Oscillator
        Method: Verlet (3rd Order)
-------------------------------------------------------------- */

using namespace std;
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>

const double dx =0.0001;                  //time step, can change to 0.001
const double tfin=50.0;                 //the maximum time we consider
const double w02=1.0;
const double xbot=-5.0;
const double xtop=5.0;


double fx(double, double, double, double);           //prototype
double Vx(double, double);                           //prototype

double abs(double x) {return sqrt(x*x);}


int main(){
    double xolder, xold, xnew, psiolder, psiold,  psinew, vold,Eincr,E,a;
    double k1y, k2y, k1v, k2v;
    Eincr = 0.001;

    for(a=0.0;a<3;a++){
        cout << setw(20) << "First three eigenvalue energies for a =" << setw(4) << a << " " << endl;

        for(E=0.0;E<4.0;E+=Eincr){
    xold=xbot;
    psiold=exp(-sqrt(2.0*(Vx(xold,a)-E))*abs(xold));
  //              ofstream prt1("qoscx4042080475.dat");

//initial conditions
vold=0.0;
xnew=xold+dx;

//RK Evolution for 1st Step
k1y=vold*dx;
k1v=fx(psiold,xold,E,a)*dx;
k2y=(vold+k1v/2.0)*dx;
k2v=fx(psiold+k1y/2.0,xold+dx/2.0,E,a)*dx;
psinew=psiold+k2y;

//Verlet Evolution
while(xnew<xtop){
xolder = xold;
psiolder = psiold;
psiold = psinew;
xold = xnew;


psinew = 2.0*psiold - psiolder + fx(psiold,xold,E,a)*dx*dx;
xnew=xold+dx;


//prt1 << setw(12) << xolder << " " << setw(12) <<  psiolder*1e19 << " "  << endl;

    }

if(abs(psiolder*2e7)<10){cout << setw(12) << setprecision(8) << E  << endl;}

/*
if(lam<2){    
if(E<1){if(abs(psiolder*1e19)<10){cout << setw(12) << setprecision(8) << E  << endl;}}
if(1<=E<=2){if(abs(psiolder*1e19)<3){cout << setw(12) << setprecision(8) << E  << endl;}}
if(2<=E<=3){if(abs(psiolder*1e19)<0.05){cout << setw(12) << setprecision(8) << E  << endl;}}}

if(lam>2){
if(E<1){if(abs(psiolder*1e19)<1){cout << setw(12) << setprecision(8) << E  << endl;}}
if(1<=E<=2){if(abs(psiolder*1e19)<0.5){cout << setw(12) << setprecision(8) << E  << endl;}}
if(2<=E<=3){if(abs(psiolder*1e19)<0.005){cout << setw(12) << setprecision(8) << E  << endl;}}}

}
if(a==2.5){a=2.0;}
*/
}
}



}

/*--------------------------------------------------------------
   Model dependent part
-------------------------------------------------------------- */
const double w0=1.0;
const double m=1.0;

double fx(double psi, double x,double E,double a){
  return 2.0*(Vx(x,a)-E)*psi;
}

double Vx(double x, double a){
 if(x>0){return 0.5*m*w0*w0*(x-a)*(x-a);}else if(x==0){return 0.5*m*w0*w0*(x-a)*(x-a);}else{ return 0.5*m*w0*w0*(x+a)*(x+a);}
}

