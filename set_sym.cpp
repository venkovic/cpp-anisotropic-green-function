#include <math.h>
#include <vector>
#include "GreenAnisotropic2D.hpp"

void medium::set_sym(int sym, std::vector<double> p) {
	isym=sym;
	if (isym==0) {
		P0=p[0];
		P1=p[1];
		R0=p[2];
		R1=p[3];
		T0=p[4];
		T1=p[5];	
		L1111=T0+2.*T1+R0*cos(4.*P0)+4.*R1*cos(2.*P1);
		L1112=R0*sin(4.*P0)+2.*R1*sin(2.*P1);
		L1122=-T0+2.*T1-R0*cos(4.*P0);
		L1212=T0-R0*cos(4.*P0);
		L2212=-R0*sin(4.*P0)+2.*R1*sin(2.*P1);
		L2222=T0+2.*T1+R0*cos(4.*P0)-4.*R1*cos(2.*P1);
	}
	else if (isym==1) {
		Th=p[0];
		K =p[1];
		R0=p[2];
		R1=p[3];
		T0=p[4];
		T1=p[5];
		L1111=T0+2*T1+pow(-1,K)*R0*cos(4*Th)+4*R1*cos(2*Th);
		L1112=-pow(-1,K)*R0*sin(4*Th)-2*R1*sin(2*Th);
		L1122=-T0+2*T1-pow(-1,K)*R0*cos(4*Th);
		L1212=T0-pow(-1,K)*R0*cos(4*Th);
		L2212=pow(-1,K)*R0*sin(4*Th)-2*R1*sin(2*Th);
		L2222=T0+2*T1+pow(-1,K)*R0*cos(4*Th)-4*R1*cos(2*Th);
	}
	else if (isym==2) {
		Th=p[0];
		R1=p[1];
		T0=p[2];
		T1=p[3];
		L1111=T0+2*T1+4*R1*cos(2*Th);
		L1112=-2*R1*sin(2*Th);
		L1122=-T0+2*T1;
		L1212=T0;
		L2212=-2*R1*sin(2*Th);
		L2222=T0+2*T1-4*R1*cos(2*Th);
	}
	else if (isym==3) {
		Th=p[0];
		R0=p[1];
		T0=p[2];
		T1=p[3];
		L1111=T0+2*T1+R0*cos(4*Th);
		L1112=-R0*sin(4*Th);
		L1122=-T0+2*T1-R0*cos(4*Th);
		L1212=T0-R0*cos(4*Th);
		L2212=R0*sin(4*Th);
		L2222=T0+2*T1+R0*cos(4*Th);	
	}
	else if (isym==4) {
		T0=p[0];
		T1=p[1];
		kappa=2.*T1;
		mu=T0;
		L1111=kappa+mu;
		L1112=0.;
		L1122=kappa-mu;
		L1212=mu;
		L2212=0.;
		L2222=kappa+mu;
	}
	else if (isym==5) {
		kappa=p[0];
		mu=p[1];
		L1111=kappa+mu;
		L1112=0.;
		L1122=kappa-mu;
		L1212=mu;
		L2212=0.;
		L2222=kappa+mu;
	}
	L2211=L1122;
	L1121=L1112;L1211=L1112;L2111=L1112;
	L2221=L2212;L1222=L2212;L2122=L2212;
	L1221=L1212;L2112=L1212;L2121=L1212;
}

double medium::dkN1(double t, int k, int i, int j) {
	if (isym==0) {
		return dkN1_anisotropic(t,k,i,j);
	}
	else if (isym==1) {
		return dkN1_orthotropic(t,k,i,j);
	}
	else if (isym==2) {
		return dkN1_r0_orthotropic(t,k,i,j);
	}
	else if (isym==3) {
		return dkN1_square_symmetric(t,k,i,j);
	}	
	else if ((isym==4)||(isym==5)) {
		return dkN1_isotropic(t,k,i,j);
	}
	return 0;
}

double medium::dkN2(double t, int k, int i, int j) {
	if (isym==0) {
		return dkN2_anisotropic(t,k,i,j);
	}
	else if (isym==1) {
		return dkN2_orthotropic(t,k,i,j);
	}
	else if (isym==2) {
		return dkN2_r0_orthotropic(t,k,i,j);
	}
	else if (isym==3) {
		return dkN2_square_symmetric(t,k,i,j);
	}	
	else if ((isym==4)||(isym==5)) {
		return dkN2_isotropic(t,k,i,j);
	}
	return 0;
}
