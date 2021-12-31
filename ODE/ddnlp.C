/*--------------------------------------------------------------
 ODE: Solution
        Model:  Damped Driven Non-Linear Pendulumn
        Method: Verlet (3rd Order)
-------------------------------------------------------------- */

using namespace std;
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>

const double dt =0.05;                  //time step
const double tfin=50.0;                 //the maximum time we consider
const double OM = 1.00;
const double a0=1.00;
const double gam=1.00;
const double w02=1.0;
const double PI=acos(-1.0);
const double m=1.00;

ofstream prt("DDNLP02.dat");

double ft(double y, double v, double t) { return -w02*w02*sin(y)-gam*v+a0*sin(OM*t);}

int main(){
        double told, tolder, yold, yolder, vold, volder, Eold, Eolder, tnew, ynew, vnew, Enew;
        double k1y, k1v, k2y, k2v;

//initial conditions
	told=0.0; yold=0.20; vold=0.0;
//	told=0.0; yold=1.00; vold=0.0;


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

	prt << setw(12) << told << " " << setw(12) << yold << " " << setw(12)
	<< vold << " " << endl;

	}
}

