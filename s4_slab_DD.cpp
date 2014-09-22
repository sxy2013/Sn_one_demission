//============================================================================
// Name        : s4_slab_DD.cpp
// Author      : kevin
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cmath>
#include <omp.h>
using namespace std;

// assign angles and weights
double mu[]={-0.8611363115,-0.3399810435,0.3399810435,0.8611363115};
double wt[]={0.3478548451,0.6521451549,0.6521451549,0.3478548451};

int cal_psi(double EDGE[],int n_EDGE,double NFM[],int n_NFM,double SigT[],double SigS[],double RegMat[],double Source[])
{
	int totNFM=0;
	for(int i=0;i<n_NFM;i++){
		totNFM+=NFM[i];
	}
	double psi[totNFM+1][4];
	double phi[totNFM];
	for(int i=0;i<totNFM;i++)
		phi[i]=0.0;
	double S[totNFM][4];
	double Q[totNFM][4];
	int fmmid[totNFM];
	double Delta[totNFM];
	double A[totNFM][4];
	double B[totNFM][4];

	//compute discretization

	int sum=0;
	for(int i=0;i<n_NFM;i++){
		for(int k=sum;k<sum+NFM[i];k++){
			Delta[k]=(EDGE[i+1]-EDGE[i])/NFM[i];
			fmmid[k]=RegMat[i];
			for (int j=0;j<4;j++)
				S[k][j]=Source[i];
		}
		sum=sum+NFM[i];
	}


	//precompute coefficients, following Eqs. 20.16-20.20

	int alpha=0;
	for (int i=0;i<totNFM;i++){
		int m=fmmid[i];
		for (int n=0;n<4;n++){
			int smu;
		    if(mu[n]>0)
		    	smu=1;
			else
				smu=-1;

		    double demon=2*mu[n]+smu*(1+smu*alpha)*SigT[m]*Delta[i];

		    A[i][n]=(2*mu[n]-smu*(1-smu*alpha)*SigT[m]*Delta[i])/demon;

		    B[i][n]=smu*2*Delta[i]/demon;
			}
		}


	//convergence parameters

	double eps_phi=pow(10,-25);
	int max_it=300;
	double err_phi=1.0;
	int it=0;

	//begin source iterations

	while(err_phi>eps_phi&&it<=max_it){
		double phi0[totNFM];
		for (int i=0;i<totNFM;i++){
			phi0[i]=phi[i];
		}
		for (int i=0;i<totNFM;i++){
			for (int j=0;j<4;j++){
				Q[i][j]=S[i][j]+0.5*SigS[fmmid[i]]*phi[i];
			}
		}
		// Perform sweeps
		// right to left
		for (int i=totNFM-1;i>=0;i--){
			for (int j=0;j<2;j++){
				psi[i][j]=A[i][j]*psi[i+1][j]+B[i][j]*Q[i][j];
			}
		}
		//left to right
		for (int i=0;i<totNFM;i++){
			for(int j=3;j<=4;j++){
				psi[i+1][j]=A[i][j]*psi[i][j]+B[i][j]*Q[i][j];
			}
		}
		// ---Insert Acceleration Here---
		// Update phi (need cell-centered psi, so use Eq.20.15)
		for (int i=0;i<totNFM;i++){
			for (int j=0;j<4;j++){
				double sum=0.0;
				sum=sum+wt[j]*0.5*(psi[i][j]+psi[i+1][j]);
				phi[i]=sum;
			}
		}
		double tmp[totNFM];
		for (int i=0;i<totNFM;i++){
			tmp[i]=abs(phi[i]-phi0[i])/phi0[i];
		}
		double max=tmp[0];
		for (int i=1;i<totNFM;i++){
			if (tmp[i]>max)
				max=tmp[i];
		}
		err_phi=max;
		it=it+1;
	}


	if (it<max_it)
		cout<<"Coveraged in "<<it<<" iterations"<<endl;
	else
		cout<<"Failed to coverage."<<endl;
	return 0;
}

int main(){
	double edge[]={0,2,4,6,9,10,16};
	double nfm[]={100,361,150,300,264};
	double regmat[]={0,1,2,3,4,4};
	double sigt[]={0.1,0.5,0.3,0.45,0.31};
	double sigs[]={0.05,0.2,0.11,0.32,0.22};
	double source[]={10,52,36,31,46,11};
	/*
    for (int i=0;i<101;i++){
    	edge[i]=i;
    }
	for (int i=0;i<100;i++){
		nfm[i]=50+100*i;
		regmat[i]=i;
		sigt[i]=0.1+i/20;
		sigs[i]=0.5*(0.1+i/20);
		source[i]=i%4+i*20;
	}
*/
	cal_psi(edge,7,nfm,6,sigt,sigs,regmat,source);
	return 0;
}
