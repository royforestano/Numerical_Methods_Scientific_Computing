//---------------------
using namespace std;
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <fstream>

ofstream prt("p1.dat");

const double r=0.5; //upper limit is r<0.675
const int N=1000;
const int L=40;

double g(double x) { return r*sin(3.14159*x);}

int main(){

double xn;
int i;
xn=0.9;

//cout << xn << " " << 0.0 << endl;
for(i=0;i<L;i++){
prt << xn << " " << g(xn)
<< endl;
xn=g(xn);
prt << xn << " " << xn << endl;
}
}


