//---------------------
// This will print the behavior of r with respect to x
// for the sine map x_n+1 = rsin(pi*x_n)
//---------------------
using namespace std;
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <fstream>>

ofstream prt("problem1.dat");

const int N=1000;
const int L=64;
const double rinc =0.0001;

double g(double x,double r) { return r*sin(3.14159*x);}

int main(){

double x,r;
int i;
x=0;

r=0.7;
for(r=0.7;r<=0.9;r+=rinc){
	for(i=0;i<N;i++){x=g(x,r);}
	for(i=0;i<L;i++){
	prt << r << " " << x << endl;
	x=g(x,r);
	prt << r  << " " << x << endl;
	x=x+0.00005;
	}
}
}




