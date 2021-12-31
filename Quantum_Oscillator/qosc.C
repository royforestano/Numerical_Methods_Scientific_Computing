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


double fx(double, double, double);           //prototype
double Vx(double);                           //prototype

double abs(double x) {return sqrt(x*x);}


int main(){
    double xolder, xold, xnew, psiolder, psiold,  psinew, vold,Eincr,E;
    double k1y, k2y, k1v, k2v;

//Energies for a = 0.0
//	E=0.5;		//modify psiolder by 1, psioldersq by 2
	E=1.5;		//modify psiolder by 1.5, psioldersq by 3

//Energies for a = 1.0
//	E=1.5;		//modify psiolder by 0.3, psioldersq by 0.5
//	E=2.197;		//modify psiolder by 0.75, psioldersq by 1.0

//Energies for a = 2.0
//	E=0.518;            //modify psiolder by 0.01, psioldersq by 0.005
//      E=1.367;            //modify psiolder by 0.005, psioldersq by 0.005

    xold=xbot;
    psiold=exp(-sqrt(2.0*(Vx(xold)-E))*abs(xold));

//	Changes according to what value of E we consdier
		ofstream prt1("qosca015.dat");
		ofstream prt2("qosca015sq.dat");


//initial conditions
vold=0.0;
xnew=xold+dx;

//RK Evolution for 1st Step
k1y=vold*dx;
k1v=fx(psiold,xold,E)*dx;
k2y=(vold+k1v/2.0)*dx;
k2v=fx(psiold+k1y/2.0,xold+dx/2.0,E)*dx;
psinew=psiold+k2y;

//Verlet Evolution
while(xnew<xtop){
xolder = xold;
psiolder = psiold;
psiold = psinew;
xold = xnew;


psinew = 2.0*psiold - psiolder + fx(psiold,xold,E)*dx*dx;
xnew=xold+dx;


	prt1 << setw(12) << xolder << " " << setw(12) <<  psiolder*1.5*1e5 + E << " "  << endl;
	prt2 << setw(12) << xolder << " " << setw(12) <<  (psiolder*3.0*1e5)*(psiolder*3.0*1e5) + E << " "  << endl;

    }

}



/*--------------------------------------------------------------
   Model dependent part
-------------------------------------------------------------- */
const double w0=1.0;
const double m=1.0;
const double a=0.0;
//const double a=1.0;
//const double a=2.5;

double fx(double psi, double x,double E){
  return 2.0*(Vx(x)-E)*psi;
}

double Vx(double x){
 if(x>0){return 0.5*m*w0*w0*(x-a)*(x-a);}else if(x==0){return 0.5*m*w0*w0*(x-a)*(x-a);}else{ return 0.5*m*w0*w0*(x+a)*(x+a);}
}

