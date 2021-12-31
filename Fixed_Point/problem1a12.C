//-----------------------------------
// This program will use the Newton-Raphson method to find the stable solutionsa
// of tanh(ax) - x = 0 for values of a=1.2. I hope to first generate the sequence
// to find when g(x) - x = 0. Then, take this value of x and evaluate f(x) = g'(x) = 'alpha'
// This should  retun the value 'alpha' for the approximate derivative at the fixed point.
// If it is less than one, the fixed point is stable. If the value is greater than one,
// the fixed point is unstable.
//-----------------------------------------------------------------
using namespace std;
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>

const double a=1.2; //

int N=10;

double f(double x) { return tanh(a*x) - x; }

double g(double x) { return (a*x*(1/cosh(a*x)/cosh(a*x)) - tanh(a*x))/(a*(1/cosh(a*x)/cosh(a*x))-1); }

double abs(double x){ return sqrt(x*x);}

int main() {
double x;
int i,k;

k=0;
while(k<4){

if(k==0){ x=-2;}
if(k==1){ x=-2;}
if(k==2){ x=0.2;}
if(k==3){ x=0.3;}

cout << " " << endl;
if(abs(x)>1){ cout << "     " << "We will now find a stable fixed point by starting with x=" << x << ".\n";}
if(abs(x)<1){ cout << "     " << "We will now find an unstable(runs away at small x)/" << endl << "     " << "stable(runs toward at super small x) fixed point at x=0 by starting with x=" << x << ".\n";}

for(i=0;i<N;i++){ 	cout << setw(6) << i << setw(22) << setprecision(12) << x
			<< setw(22) << f(x) << setw(22) <<  g(x)-x << endl; x=g(x);}

cout << "     " << "There is a stable fixed point at x = " << x << " for the function f(x) = tanh(1.2*x) -x.\n";
cout << " " << endl;

k++;
}

}

