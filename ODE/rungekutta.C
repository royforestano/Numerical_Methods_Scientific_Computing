/*--------------------------------------------------------------
 ODE: Solution
        Model:  Simple Harmonic Oscillator
        Method: Runge-Kutta (2nd Order)
-------------------------------------------------------------- */

using namespace std;
#include <iostream>
#include <iomanip>

const double om = 1.00;
const double dt = 0.05;                  //time step
const double tfin=50.0;                 //the maximum time we consider
const double w02=1.0;


double ft(double,double,double);        //prototype
double Et(double,double);               //prototype

int main(){
	double told, yold, vold, Eold, tnew, ynew, vnew, Enew;
	double k1y, k1v, k2y, k2v;

//initial conditions
told = 0.0;
yold=1.0;
vold=0.0;
cout << setw(12) << told << " " << setw(12) << yold << " " << setw(12)
     << vold << " " << setw(12) << Et(yold,vold)  << endl;

//RK evolution
k1y=vold*dt;
k1v=ft(yold,vold,told)*dt;
k2y=(vold+k1v/2.0)*dt;
k2v=ft(yold+k1y/2.0, vold+k1v/2.0,told*dt/2.0)*dt;
ynew = yold+k2y;
vnew=vold+k2v;
tnew=told+dt;

while(tnew<tfin){
told=tnew;
vold=vnew;
yold=ynew;

cout << setw(12) << told << " " << setw(12) << yold << " " << setw(12) 
<< vold << " " << setw(12) << Et(yold,vold)  << endl;

k1y=vold*dt;
k1v=ft(yold,vold,told)*dt;
k2y=(vold+k1v/2.0)*dt;
k2v=ft(yold+k1y/2.0, vold+k1v/2.0,told*dt/2.0)*dt;
ynew = yold + k2y;
vnew = vold + k2v;
tnew=told+dt;
Enew=Et(ynew,vnew);

}

}


/*--------------------------------------------------------------
   Model dependent part
-------------------------------------------------------------- */
const double w0=1.0;
const double m=1.0;

double ft(double x, double v, double t){
  return -w0*w0*x  ;
}
double Et(double x, double v){
  return 0.5*m*(v*v+w02*x*x);
}




