//-------------------------------------------------
// Solution of Symmetric 2x2 Eigenvalue problem 
// using Jacobi Method
//-------------------------------------------------
using namespace std;
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>

int main(){
 double a,b,d,t1,t2,theta,t,c,s,eig1,eig2;

// ifstream read("datain.dat");
// read >> a >> b >> d;

a=1;
b=2;
d=1;

 theta=(d-a)/(2.0*b); 
 t1=-theta+sqrt(theta*theta+1.0);
 t2=-theta-sqrt(theta*theta+1.0);

 if(fabs(t1) < fabs(t2)){
     t=t1;}
  else{
     t=t2;}

 c=1.0/sqrt(1.0+t*t);
 s=c*t;
 eig1=c*c*a-2.0*c*s*b+s*s*d;
 eig2=c*c*d+2.0*c*s*b+s*s*a;

//ofstream prt("dataout.dat");
 
cout << setw(16) <<  "eig1=" <<  eig1 <<  setw(16) << "eig2=" << eig2 << endl; 
cout << setw(16) << "vec1=" <<  "[" << c << ", " << -s << "]" << endl;
cout << setw(16) << "vec2=" <<  "[" << s << ", " <<  c << "]" << endl;  

}


