//---------------------
using namespace std;
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <fstream>

ofstream prt("nrm.dat");

const double a=1.2;
const int N=10;

double g(double x) { return (a*x*(1/cosh(a*x)/cosh(a*x))-tanh(a*x))/(a*(1/cosh(a*x)/cosh(a*x))-1.0); }

int main(){ 

double xn;
int i;

xn=0.3;
//cout << xn << " " << 0.0 << endl;
for(i=0;i<N;i++){
prt << xn << " " << g(xn) 
<< endl;
xn=g(xn);
prt << xn << " " << xn << endl;
}
}
