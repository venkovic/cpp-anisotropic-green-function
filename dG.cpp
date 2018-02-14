#include <math.h>
#include <vector>
#include "GreenAnisotropic2D.hpp"

using namespace std;

double medium::dh(int n, int k, vector<int> i, double th) {
	if (n==1) {
		return dkh1i(k,i,th);
	}
	else if (n>1) {
		double dkhni=0.;		
		int ival=i.back();
		i.pop_back();
		for (int s=0;s<k+1;++s) {
			dkhni+=Binom(k,s)*((n-1.)*dh(n-1,k-s,i,th)*dknvec(s,ival,th)-dh(n-1,k-s+1,i,th)*dknvec(s+1,ival,th));
		}		
		return dkhni;
	}
	return 0;
}
	
double medium::h(int n, vector<int> i, double th) {
	vector<double> nvec(2), mvec(2);	
	nvec={cos(th),sin(th)};
	mvec={-sin(th),cos(th)};
	if (n==1) {
		return dkh1i(0,i,th);
	}
	else if (n>1) {
		int ival=i.back();
		i.pop_back();
		double hni=(n-1)*h(n-1,i,th)*nvec[ival-1]-dh(n-1,1,i,th)*mvec[ival-1];
		return hni;
	}
	return 0;
}

double medium::dkh1i(int k, vector<int> i, double th) {
	vector<double> nvec(2), mvec(2);
	nvec={cos(th),sin(th)};
	mvec={-sin(th),cos(th)};
	int ival=i.back();
	i.pop_back();	
	if (k==0) {
		double h1i=H[i[0]-1][i[1]-1]*nvec[ival-1];
		h1i+=(dkN1(th,0,i[0],1)*H[1-1][i[1]-1]+dkN1(th,0,i[0],2)*H[2-1][i[1]-1]+dkN2(th,0,i[0],1)*S[i[1]-1][1-1]+dkN2(th,0,i[0],2)*S[i[1]-1][2-1])*mvec[ival-1];	
		return h1i;
	}
	else if (k>0) {
		double dkh1i=H[i[0]-1][i[1]-1]*dknvec(k,ival,th);
		for (int s=0;s<k+1;++s) {
			dkh1i+=Binom(k,s)*(dkN1(th,k-s,i[0],1)*H[1-1][i[1]-1]+dkN1(th,k-s,i[0],2)*H[2-1][i[1]-1]+dkN2(th,k-s,i[0],1)*S[i[1]-1][1-1]+dkN2(th,k-s,i[0],2)*S[i[1]-1][2-1])*dkmvec(s,ival,th);
		}
		return dkh1i;
	}
	return 0;
}

double medium::dknvec(int k, int i, double th) {
	vector<double> nvec(2), mvec(2);
	nvec={cos(th),sin(th)};
	mvec={-sin(th),cos(th)};
	return pow(-1,k*(k+3)/2)*pow(nvec[i-1],(pow(-1,k)+1)/2)*pow(mvec[i-1],(pow(-1,k+1)+1)/2);
}

double medium::dkmvec(int k, int i, double th) {
	vector<double> nvec(2), mvec(2);
	nvec={cos(th),sin(th)};
	mvec={-sin(th),cos(th)};
	return pow(-1,(1./2.)*(k+1)*(k+4))*pow(mvec[i-1],(1./2.)*pow(-1,k+2)+1./2.)*pow(nvec[i-1],(1./2.)*pow(-1,k+1)+1./2.);	
}

double medium::h_DP(int n, vector<int> i, double th) {
	vector<double> nvec(2), mvec(2);	
	nvec={cos(th),sin(th)};
	mvec={-sin(th),cos(th)};
	//
	vector<double> d0hk(n);
	for (int k=1;k<n+1;++k) {
		for (int rr=0;rr<n-k+1;++rr) {
			int r=n-k-rr;
			for (int s=0;s<r+1;++s) {
				if (s==0) {
					if (k==1) {
						d0hk[r+k-1]=H[i[0]-1][i[1]-1]*dknvec(r,i[1+k],th)\
						            +(dkN1(th,r,i[0],1)*H[1-1][i[1]-1]+dkN1(th,r,i[0],2)*H[2-1][i[1]-1]+dkN2(th,r,i[0],1)*S[i[1]-1][1-1]+dkN2(th,r,i[0],2)*S[i[1]-1][2-1])*dkmvec(0,i[1+k],th);
					}
					else {
						d0hk[r+k-1]=(k-1.)*d0hk[r+k-2]*dknvec(0,i[1+k],th)-d0hk[r+k-1]*dknvec(1,i[1+k],th);
					}
				}
				else {
					if (k==1) {
						d0hk[r+k-1]+=Binom(r,s)*(dkN1(th,r-s,i[0],1)*H[1-1][i[1]-1]+dkN1(th,r-s,i[0],2)*H[2-1][i[1]-1]+dkN2(th,r-s,i[0],1)*S[i[1]-1][1-1]+dkN2(th,r-s,i[0],2)*S[i[1]-1][2-1])*dkmvec(s,i[1+k],th);
					}
					else {
						d0hk[r+k-1]+=Binom(r,s)*((k-1.)*d0hk[r-s+k-2]*dknvec(s,i[1+k],th)-d0hk[r-s+k-1]*dknvec(s+1,i[1+k],th));
					}
				}				
			}		
		}
		// At this stage, r in [0,n-k] => d0hk[r+k-1] = d^{(r)}_{\theta}[h^{k}_{ijk_1\dots k_k}(\theta)]	
	}
	// At this stage, k in [1,n] => d0hk[k-1] = h^{k}_{ijk_1\dots k_k}(\theta)
	return d0hk[n-1];

}

double medium::dnGi(int n, vector<int> i, double r, double th) {
	double gn=pow(-r,-n)/M_PI;
	vector<double> nvec(2), mvec(2);
	nvec={cos(th),sin(th)};
	mvec={-sin(th),cos(th)};
	if (n>=1) {
		// 
		// Recursive implementation:
		//return gn*h(n,i,th)/2.;
		//
		// Bottom-up dynamic progamming implementation:
		return gn*h_DP(n,i,th)/2.;
	}	
	return 0;
}

