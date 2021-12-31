
//-----------------------------------
// This program will use the traditional method to find the stable solutions a
// of tanh(ax) = x for values of a=0.8 and a=1.2. I hope to first generate the sequence
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
const double b=1.2;

double f(double x) { return tanh(a*x); }
double g(double x) { return tanh(b*x); }

int main() {

double x;
int i;

x=0.001;

cout << endl;
cout << endl;
cout<< "     " << "This program will produce the fixed point resullts for the traditional iterative method for a=08 and a=1.2 for the logistic map." << endl;
cout << endl;

for(i=0;i<100;i++) {

cout << setw(6) << i << setw(22) << setprecision(12) << x
<< setw(22) << f(x) << setw(22)  << endl; x=f(x);

  }
if(x<1e-10){x=0;}
cout << "There is a stable point at x = " << x
<< " for the function f(x) = tanh(0.8*x) -x.\n";

cout << "We will now show the solutions for a=1.2 in the traditional method.\n";


x=0.55;

for(i=0;i<100;i++) {

cout << setw(6) << i << setw(22) << setprecision(12) << x
<< setw(22) << g(x) << setw(22)  << endl; x=g(x);

  }

cout << "There is a stable point at x = " << x
<< " for the function f(x) = tanh(1.2*x) -x.\n";



x=-0.55;

for(i=0;i<100;i++) {

cout << setw(6) << i << setw(22) << setprecision(12) << x
<< setw(22) << g(x) << setw(22)  << endl; x=g(x);

  }

cout << "There is a stable point at x = " << x
<< " for the function f(x) = tanh(1.2*x) -x.\n";

cout << "Starting from a small positive number we find.\n";

x=0.01;

for(i=0;i<20;i++) {

cout << setw(6) << i << setw(22) << setprecision(12) << x
<< setw(22) << g(x) << setw(22)  << endl; x=g(x);

  }

cout << "Starting from a small negatuve number we find.\n";

x=-0.01;

for(i=0;i<20;i++) {

cout << setw(6) << i << setw(22) << setprecision(12) << x
<< setw(22) << g(x) << setw(22)  << endl; x=g(x);

  }


cout << "There is an unstable point at x = " << "0" << " for the function f(x) = tanh(1.2*x) -x.\n";
cout << "One can see when we start significantly close to 0, "
<< "the system trends away from x=0 to the stable fixed points found above.\n";


}

