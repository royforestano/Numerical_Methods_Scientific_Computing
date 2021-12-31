/*--------------------------------------------------------------
 ODE: Solution
        Model:  Mackey Glass
        Method: Runge-Kutta (4th Order)
-------------------------------------------------------------- */

using namespace std;
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>


const double beta = 0.20;
const double gam = 0.10;
const double tau = 17.0;
const int    nmg = 10;
const double dt = 0.10;                  //time step
const double tfin=1000.0;                 //the maximum time we consider
const int    ntmax=int(tfin/dt)+1;
const int    ntau =int(tau/dt);
double u[ntmax];

ofstream prt1("mglass.dat");


double ft(double,double,double);        //prototype

int main(){
	double told, tnew, yold, ynew;
	double k1, k2, k3, k4;
	int nt,i;

//initial conditions
told = 0.0;
yold=0.8;
nt = 0;
u[nt]=yold;


//RK evolution
k1 = ft(told,yold,0.0);
k2=ft(told+dt/2.0,yold+dt*k1/2.0,0.0);
k3=ft(told+dt/2.0,yold+dt*k1/2.0,0.0);
k4 = ft(told+dt,yold+dt*k3,0.0);
ynew = yold+dt*(k1+2.0*k2+2.0*k3+k4)/6.0;
tnew=told+dt;
nt++;
u[nt]=ynew;

while(tnew<tau){
told=tnew;
yold=ynew;

prt1 << setw(12) << told << " " << setw(12) << yold << " " << setw(12) << 0.0 << endl; 

k1 = ft(told,yold,0.0);
k2=ft(told+dt/2.0,yold+dt*k1/2.0,0.0);
k3=ft(told+dt/2.0,yold+dt*k1/2.0,0.0);
k4 = ft(told+dt,yold+dt*k3,0.0);
ynew = yold+dt*(k1+2.0*k2+2.0*k3+k4)/6.0;
tnew=told+dt;
nt++;
u[nt]=ynew;

}


while(tnew<tfin){
told=tnew;
yold=ynew;

prt1 << setw(12) << told << " " << setw(12) << yold << " " << setw(12) << u[nt-ntau] << endl;

k1 = ft(told,yold,u[nt-ntau]);
k2=ft(told+dt/2.0,yold+dt*k1/2.0,u[nt-ntau]);
k3=ft(told+dt/2.0,yold+dt*k1/2.0,u[nt-ntau]);
k4 = ft(told+dt,yold+dt*k3,u[nt-ntau]);
ynew = yold+dt*(k1+2.0*k2+2.0*k3+k4)/6.0;
tnew=told+dt;
nt++;
u[nt]=ynew;

}

}


double ft(double t, double u, double up){
  return beta*up/(1.0+exp(nmg*log(up)))-gam*u;
}


