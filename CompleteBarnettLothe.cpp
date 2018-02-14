#include <math.h>
#include <vector>
#include "GreenAnisotropic2D.hpp"

using namespace std;

double medium::trapzd(double a, double b, int n, int Ni, int I, int J) {
	double x,tnm,sum,del;
	sum=0.;
	static double s;
	int it,j;
	
	if (n==1) {
		if (Ni==1) {
			return (s=.5*(b-a)*(dkN1(a,0,I,J)+dkN1(b,0,I,J)));
		}
		else if (Ni==2) {
			return (s=.5*(b-a)*(dkN2(a,0,I,J)+dkN2(b,0,I,J)));
		}
	}
	else {
		for (it=1,j=1;j<n-1;j++) it<<=1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+.5*del;
		if (Ni==1) {
			for (sum=0.,j=0;j<it;j++,x+=del) {
				sum+=dkN1(x,0,I,J);
			}
		}
		else if (Ni==2) {
			for (sum=0.,j=0;j<it;j++,x+=del) {
				sum+=dkN2(x,0,I,J);
			}
		}
		s=.5*(s+(b-a)*sum/tnm);
		return s;
	}
	return 0;
}

double medium::qtrap(double a, double b, int Ni, int I, int J) {
	const int JMAX=20;
	const double EPS=1.e-10;
	int j;
	double s,olds=0.;
	
	for (j=0;j<JMAX;j++) {
		s=trapzd(a,b,j+1,Ni,I,J);
		if (j>5) {
			if ((fabs(s-olds)<EPS*fabs(olds))||(s==0.&&olds==0.)) {
				return s;
			}
		}
		olds=s;
	}
	return 0;
}

void medium::get_S() {
	S={{0,0},{0,0}};
	S[0][0]=qtrap(0.,M_PI,1,1,1)/M_PI;
	S[1][1]=qtrap(0.,M_PI,1,2,2)/M_PI;
	S[0][1]=qtrap(0.,M_PI,1,1,2)/M_PI;
	S[1][0]=qtrap(0.,M_PI,1,2,1)/M_PI;
}

void medium::get_H() {
	H={{0,0},{0,0}};
	H[0][0]=qtrap(0.,M_PI,2,1,1)/M_PI;
	H[1][1]=qtrap(0.,M_PI,2,2,2)/M_PI;
	H[0][1]=qtrap(0.,M_PI,2,1,2)/M_PI;
	H[1][0]=qtrap(0.,M_PI,2,2,1)/M_PI;
}
