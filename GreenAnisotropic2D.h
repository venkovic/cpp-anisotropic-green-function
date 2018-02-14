#include <pybind11/pybind11.h>
#include <vector>	

namespace py = pybind11;

#ifndef MEDIUM_HEADER
#define MEDIUM_HEADER

struct medium{
	// Symmetry code
	int isym;
	//	isym=0 -> Anistropic medium 		// (P0,P1,R0,R1,T0,T1)
	//	isym=1 -> Orthotropic medium 		// (Th,K,R0,R1,T0,T1)
	//	isym=2 -> R0-Orthotropic medium 	// (Th,R1,T0,T1)
	//	isym=3 -> Square symmetric medium	// (Th,R0,T0,T1)
	//	isym=4 -> Polar isotropic medium 	// (T0,T1)
	//	isym=5 -> Isotropic medium 			// (m,k)
	//
	// Polar invariants
	double P0,P1,R0,R1,T0,T1;
	int K;
	//
	// Isotropic material constants 
	double kappa, mu;
	//
	// Phase angle 
	double Th;
	//
	// Stiffness components
	double L1111,L1122,L1112;
	double L2211,L2222,L2212;
	double L1211,L1222,L1212;
	double L1121;
	double L2221;
	double L2111,L2122,L2112,L1221,L2121;
	//
	// First complete Barnett-Lothe tensor integral
	std::vector<std::vector<double>> S;
	//
	// Second complete Barnett-Lothe tensor integral
	std::vector<std::vector<double>> H;
	
	
	// Sets up material symmetry
	//void set_sym(int sym, py::list p);
	void set_sym(int sym, std::vector<double> p);
	//
	// Computes first complete Barnett-Lothe tensor integral
	void get_S();
	//
	// Computes second complete Barnett-Lothe tensor integral
	void get_H();
	//
	// Computes k-th derivative of first Barnett-Lothe integrand
	double dkN1_anisotropic(double t, int k, int i, int j);
	double dkN1_orthotropic(double t,int k, int i, int j);
	double dkN1_r0_orthotropic(double t, int k, int i, int j);
	double dkN1_square_symmetric(double t, int k, int i, int j);
	double dkN1_isotropic(double t, int k, int i, int j);
	double dkN1(double t, int k, int i, int j);
	//double dkN1_isotropic();
	//
	// Computes k-th derivative of second Barnett-Lothe integrand
	double dkN2_anisotropic(double t, int k, int i, int j);
	double dkN2_orthotropic(double t,int k, int i, int j);
	double dkN2_r0_orthotropic(double t, int k, int i, int j);
	double dkN2_square_symmetric(double t, int k, int i, int j);
	double dkN2_isotropic(double t, int k, int i, int j);
	double dkN2(double t, int k, int i, int j);	
	//double dkN1_isotropic();
	//
	// Computes k-th derivative of angle dependent part of 
	// n-th gradient of Green's function
	double dh(int n, int k, std::vector<int> i, double th);
	//
	// Computes angle dependent part of k-th gradient of Green's function
	double h(int n, std::vector<int> i, double th);
	double h_DP(int n, std::vector<int> i, double th);
	//
	// Computes k-th derivative of angle dependent part of 1-st gradient the Green's function
	double dkh1i(int k, std::vector<int> i, double th);
	//
	// Computes k-th derivative of Cartesian component i of unit vector nvec
	double dknvec(int k, int i, double th);
	//
	// Computes k-th derivative of Cartesian component i of unit vector mvec, 
	// i.e. active counter-clowise Pi/2 rotation of nvec
	double dkmvec(int k, int i, double th);
	//
	// Computes Cartesian component with indices i of n-th gradient Green's function
	double dnGi(int n, std::vector<int> i, double r, double th);
	//
	// Numerical integration
	double trapzd(double a, double b, int n, int Ni, int I, int J);
	double qtrap(double a, double b, int Ni, int I, int J);
};
//
// Binomial coefficient
int Binom(int n, int k);

#endif