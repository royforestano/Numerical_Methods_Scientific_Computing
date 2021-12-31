/*--------------------------------------------------------------
 ODE: Solution 
        Model:  Simple Harmonic Oscillator 
        Method: Euler 
-------------------------------------------------------------- */

using namespace std;
#include <iostream>
#include <iomanip>


const double dt =0.10;			//time step
const double tfin=50.0;			//the maximum time we consider
const double w02=1.0;


double ft(double,double,double);	//prototype
double Et(double,double);		//prototype

int main(){
  double told, xold, vold, Eold, tnew, xnew, vnew, Enew;


//initial conditions
tnew = 0.0;
xnew=1.0;
vnew=0.0;
Enew=Et(xnew,vnew);


//Euler evolution



while(tnew<tfin){
xold=xnew;
vold=vnew;
Eold=Enew;
told=tnew;
cout << setw(10) << told << " " << setw(10) << xold << " " << setw(10) << vold << " "
     << setw(10) << Eold  << endl;

vnew=vold+dt*ft(xold,vold,told);
xnew=xold+dt*vold;
tnew=told+dt;
Enew=Et(xnew,vnew);

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

