/*--------------------------------------------------------------
 ODE: Solution
        Model:  Damped Driven Non-Linear Pendulumn (a0 Higher)
        Method: Verlet (3rd Order)
-------------------------------------------------------------- */

using namespace std;
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

			                   //time step
const double tfin=10000.0;                 //the maximum time we consider
const double OM = 2.00/3.00;
const double a0=1.20;
const double gam=0.50;
const double w02=1.0;
const double PI=acos(-1.0);
const double m=1.00;
const double TD=2.00*PI/OM;
const int    nn=100;
const double dt=TD/nn;


double ft(double y, double v, double t) { return -w02*w02*sin(y)-gam*v+a0*sin(OM*t);}

int main(){
        double told, tolder, yold, yolder, vold, volder, Eold, Eolder, tnew, ynew, vnew, Enew;
        double k1y, k1v, k2y, k2v;
	int nt;
	ofstream prt1("ddnlp0512.dat");
	ofstream prt2("ddnlp0512poinc.dat");
	ofstream prt3("ddnlp0512poinc2.dat");

	
//initial conditions
told=0.0; yold=0.20; vold=0.0;


//RK Evolution for 1st Step
k1y=vold*dt;
k1v=ft(yold,vold,told)*dt;
k2y=(vold+k1v/2.0)*dt;
k2v=ft(yold+k1y/2.0,vold+k1v/2.0,told+dt/2.0)*dt;
ynew=yold+k2y;
vnew=vold+k2v;
tnew=told+dt;

//Verlet Evolution
while(tnew<tfin){
tolder = told;
yolder = yold;
volder = vold;
told = tnew; 
yold = ynew;


ynew = (2.0*yold - yolder*(1.0-gam*dt/2.0) 
	- (w02*w02*sin(yold)-a0*sin(OM*told))*dt*dt
	)/(1.0+gam*dt/2.0);
vold=(ynew-yolder)/dt/2.0;
tnew=told+dt;
nt++;

if (told<500.0){prt1  << setw(12) << tolder << " " << setw(12) 
<< yolder << " " << setw(12) << volder << endl;}

if ((nt-2)%nn == 0){prt2  << setw(12) << yolder << " " << setw(12) << volder << endl;}

if ((nt-2)%nn == nn/2){prt3  << setw(12) << yolder << " " << setw(12) << volder << endl;}




}
}

