//---------------------------------------------------------------
//
//	
//
//---------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdlib.h>
using namespace std;

const int NS=100*1024;
const int NT=100*100;

double r(double x, double y){ return sqrt(x*x+y*y);}

main(){
	int i,j,im;
	double x,y,v,mean,av,avsq,var,sig;
	srand48(1234);
	
	av=0.0;
	avsq=0.0;
	for(j=1;j<NS;j++){
		im=0;
		for(i=1;i<NT;i++){
			x=drand48();
			y=drand48();
			v=r(x,y);
			if (v<1.0){im++;}
		}
		mean=4.0*im/double(NT);
		cout << mean << endl;
		av+=mean;
		avsq+=mean*mean;
	}
	av=av/NS;
	avsq=avsq/NS;	
	var=avsq-av*av;
	sig=sqrt(var);
	cout << "mean=" << av << " +- " << sig << endl;
}
