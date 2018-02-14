#include <math.h> 
#include "GreenAnisotropic2D.hpp"

double medium::dkN2_isotropic(double t, int k, int i, int j) {
	if ((k==0)&&(i==1)&&(j==1)) {
		return (2.*mu+kappa*(1.+cos(2.*t)))/(2.*mu*(kappa+mu));
	}
	else if ((k==0)&&(i==2)&&(j==2)) {
		return (2.*mu+kappa*(1.-cos(2.*t)))/(2.*mu*(kappa+mu));
	}
	else if ((k==0)&&(i==1)&&(j==2)) {
		return kappa*sin(2.*t)/(2.*mu*(kappa+mu));
	}
	else if ((k==0)&&(i==2)&&(j==1)) {
		return dkN2_isotropic(t,0,1,2);
	}
	else if ((k>0)&&(i==1)&&(j==1)) {
		return pow(2.,k-1.)*kappa/(mu*(kappa+mu))*cos(k*M_PI/2.+2.*t);
	}
	else if ((k>0)&&(i==2)&&(j==2)) {
		return -dkN1_isotropic(t,k,1,1);
	}
	else if ((k>0)&&(i==1)&&(j==2)) {
		return pow(2.,k-1.)*kappa/(mu*(kappa+mu))*sin(k*M_PI/2.+2.*t);
	}
	else if ((k>0)&&(i==2)&&(j==1)) {
		return dkN1_isotropic(t,k,1,2);
	}
	return 0;
};
