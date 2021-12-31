/*--------------------------------------------------------------
 ODE: Solution
        Model:  Simple Harmonic Oscillator
        Method: Verlet (3rd Order)
-------------------------------------------------------------- */

using namespace std;
#include <iostream>
#include <iomanip>


const double dt =0.05;                  //time step
const double tfin=50.0;                 //the maximum time we consider
const double w02=1.0;


double ft(double);   		        //prototype
double Et(double,double);               //prototype

int main(){
        double told, tolder, yold, yolder, vold, volder, Eold, Eolder, tnew, ynew, vnew, Enew;
        double k1y, k1v, k2y, k2v;

//initial conditions
told=0.0; yold=2.0; vold=0.0;


//RK Evolution for 1st Step
k1y=vold*dt;
k1v=ft(yold)*dt;
k2y=ft(vold+k1v/2.0)*dt;
k2v=ft(yold+k1y/2.0)*dt;
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
vold = vnew;

ynew = 2.0*yold - yolder + ft(yold)*dt*dt;
vold=(ynew-yolder)/dt/2.0;
tnew=told+dt;

Enew=Et(ynew,vnew);

cout << setw(12) << told << " " << setw(12) << yold << " " << setw(12)
<< vold << " " << setw(12) << Et(yold,vold)  << endl;



}

}


/*--------------------------------------------------------------
   Model dependent part
-------------------------------------------------------------- */
const double w0=1.0;
const double m=1.0;

double ft(double x){
  return -w0*w0*x  ;
}
double Et(double x, double v){
  return 0.5*m*(v*v+w02*x*x);
}


