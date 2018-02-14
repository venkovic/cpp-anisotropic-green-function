#include <math.h> 
#include "GreenAnisotropic2D.hpp"

double medium::dkN1_isotropic(double t, int k, int i, int j) {
	if ((k==0)&&(i==1)&&(j==1)) {
		return kappa*sin(2.*t)/(kappa+mu);
	}
	else if ((k==0)&&(i==2)&&(j==2)) {
		return -dkN1_isotropic(t,0,1,1);
	}
	else if ((k==0)&&(i==1)&&(j==2)) {
		return -(mu+kappa*cos(2.*t))/(kappa+mu);
	}
	else if ((k==0)&&(i==2)&&(j==1)) {
		return (mu-kappa*cos(2.*t))/(kappa+mu);
	}
	else if ((k>0)&&(i==1)&&(j==1)) {
		return pow(2.,k)*kappa/(kappa+mu)*sin(k*M_PI/2.+2.*t);
	}
	else if ((k>0)&&(i==2)&&(j==2)) {
		return -dkN1_isotropic(t,k,1,1);
	}
	else if ((k>0)&&(i==1)&&(j==2)) {
		return -pow(2.,k)*kappa/(kappa+mu)*cos(k*M_PI/2.+2.*t);
	}
	else if ((k>0)&&(i==2)&&(j==1)) {
		return dkN1_isotropic(t,k,1,2);
	}
	return 0;
};
