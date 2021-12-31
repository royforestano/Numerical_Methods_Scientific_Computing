//-----------------------------------
// This program will use the Newton-Raphson method to find the stable solutionsa
// of tanh(ax) - x = 0 for values of a=0.8. I hope to first generate the sequence
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

const double a=0.8;

int N=10;

double f(double x) { return tanh(a*x) - x; }

double g(double x) { return (a*x*(1/cosh(a*x)/cosh(a*x)) - tanh(a*x))/(a*(1/cosh(a*x)/cosh(a*x))-1); }

int main() {

double x;
int i;

//cout << "Enter an initial x value: \n";
//cin >> x;

x=8;

for(i=0;i<N;i++) {

cout << setw(6) << i << setw(22) << setprecision(12) << x 
<< setw(22) << f(x) << setw(22) <<  g(x)-x << endl; x=g(x);
 
  }

cout << "There is a stable point at x = " << x
<< " for the function f(x) = tanh(0.8*x) -x.\n";


}

