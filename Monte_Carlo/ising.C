//***********************************************************************
//
//HW6: Metropolis Monte Carlo
// Develop a Markpv Chain of configurations of Ising spins that become
// distributed according to the Boltzmann probabiliy Distribution
//
//**********************************************************************
using namespace std;
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdlib.h>

ofstream prt("isingT.dat");

const int N=32;
const int MCsteps=4000;
const int MCtransient=500;
const int J=1;

void init_microstate(int S[N][N]){
	int i,j;
	double r;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			r=drand48();
			if(r>0.5){S[i][j]=1;}
			else {S[i][j]=-1;}

	}}}

void print_microstate(int S[N][N]){
	int i,j;
        for(i=0;i<N;i++){for(j=0;j<N;j++){
		if(S[i][j]==1){cout << " + " ;}
		else{cout << "  ";}
		}
		 cout << endl;
	}}

int energy_microstate(int S[N][N]){
	int i,j,R,T,sum;
	sum=0;
	for(i=0;i<N;i++){
		T=i+1;	T=T%N;
		for(j=0;j<N;j++){
			R=j+1;	R=R%N;	sum+=S[i][j]*(S[i][R]+S[T][j]);
	}	}
	return -sum;
}


int DeltaE(int S[N][N],int i,int j){
        int R,T,L,B,sum;
        sum=0;
		B=i-1; B=(B+N)%N;
                T=i+1;  T=T%N;
			L=j-1;	L=(L+N)%N;
                        R=j+1;  R=R%N;  sum=(-S[i][j])*(S[i][R]+S[i][L]+S[T][j]+S[B][j]);
        return -2*sum*J;

}

int main(){
	int S[N][N],M,E,DE,i,j,kmc,accept;
	double kT,ratio,om,aveE,aveM,aveE2,aveM2;
	srand48(4321);
	kT=3.0;

	while(kT>0.9){	
	init_microstate(S);
	M=0;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			M+=S[i][j];
	}}

	//print_microstate(S);
	E=energy_microstate(S);
	//cout << "E = " << E << "M = " << M <<  endl;
	aveE=0; aveM=0; aveE2=0; aveM2=0;
	for(kmc=0;kmc<MCsteps+MCtransient;kmc++){	//loop over MC steps
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				DE=DeltaE(S,i,j);
				// cout <<  DE << endl;
				if(DE<0){ accept=1; }
				else{	ratio=exp(-DE/kT); 
					om=drand48();
					if(ratio>om){	accept=1;}
					else	{	accept=0;}
				}
				if(accept==1){
					S[i][j]=-S[i][j];
					E=E+DE;
					M=M+2*S[i][j];
					// cout << "accepted with r = " << ratio << " and om = " << endl;
					}else{
					//cout << "rejected with r = " << ratio << " and om = " << om << endl;
					}
						
			}} //end MC step
		if(kmc>MCtransient-1){
			aveE+=E;
			aveM+=M;
			aveE2+=E*E;
			aveM2+=M*M;
		}
} //loop over MC steps
	aveE=(aveE)/MCsteps;
	aveM=(aveM)/MCsteps;
	aveE2=(aveE2)/MCsteps;
	aveM2=(aveM2)/MCsteps;
// cout << "ave E per site = " << aveE/N/N << endl; 
// cout << "ave M per site = " << aveM/N/N << endl;
// cout << "ave susceptibility = " << (aveM2-aveM*aveM)/(kT*kT)/N/N << endl;
// cout << "ave heat capacity = " << (aveE2-aveE*aveE)/(kT*kT)/N/N << endl;
	prt << setw(6) << kT << " ";
	prt << setw(9) << aveE/N/N << " ";
	prt << setw(9) << aveM/N/N << " ";
	prt << setw(9) << (aveM2-aveM*aveM)/(kT*kT)/N/N << " ";
	prt << setw(9) << (aveE2-aveE*aveE)/(kT*kT)/N/N << endl;

kT=kT-0.1;

} // end while

}
