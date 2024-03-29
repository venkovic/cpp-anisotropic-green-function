## green-anisotropic

##### Python wrapped C++ code for bottom-up dynamic computation of high order  gradient components of 2D linear elastic Green functions in anistropic,  orthotropic, R0-orthotropic, square symmetric and isotropic media.



Author: Nicolas Venkovic.

email: venkovic@gmail.com.



#### Supporting documents:

- Main supporting document ([GreenAnisotropic2D.pdf](https://github.com/nvenkov1/green-anisotropic/blob/master/GreenAnisotropic2D.pdf)),
- Description of the bottom-up dynamic implementation ([Bottom-up_DP_algo.pdf](https://github.com/nvenkov1/green-anisotropic/blob/master/Bottom-up_DP_algo.pdf)).

#### Compiler and dependencies:

 - gcc.
 - Python 2.X.
 - pybind11.
 - NumPy, Matplotlib, SciPy.

#### Installation: 

```bash
$ python setup.py build_ext -i
```

#### Tests:

- Global equilibrium of random bodies:
  Verifies that the integral of the traction field of a body bounded by a random curve is in equilibrium with a concentrated load applied at the origin.
  Verification done for 10 instances of each material symmetry and 10 realizations of random curves.
- Numerical differentiation of Green functions:
  Verifies that the computation of specific components of several gradients of Green's function matches pre-calculated values.
  Verification done for 5 material instances of each material symmetry with gradients up to the eighth order.

```bash
$ python tests.py
```

#### Usage:

An example script to start working with the module is given. Run the examples as follows: 

```bash
$ python example.py
```

#### C++ content wrapped in Python:

 -  isym : Symmetry code
-  P0, P1, R0, R1, T0, T1, K : Polar invariants
-  kappa, mu : Isotropic elastic constants
-  Th : Phase angle
-  L1111, L1122, L1112 : Stiffness components
-  L2211, L2222, L2212 : Stiffness components
-  L1211, L1222, L1212 : Stiffness components
-  L1121, L2221, L2111 : Stiffness components
-  L2122, L2112, L1221 : Stiffness components
-  L2121 : Stiffness component
-  S : First complete Barnett-Lothe tensor integral
-  H : Second complete Barnett-Lothe tensor integral
-  set_sym() : Sets up material symmetry
-  get_S() : Computes first complete Barnett-Lothe tensor integral
-  get_H() : Computes second complete Barnett-Lothe tensor integral
-  dnGi() : Computes Cartesian component with indices i of n-th gradient Green's function

#### C++ content not wrapped in Python:

 -  dkN1_anisotropic() : Computes k-th derivative of first Barnett-Lothe integrand of an anisotropic medium
-  dkN2_anisotropic() : Computes k-th derivative of second Barnett-Lothe integrand of an anisotropic medium
-  dkN1_orthotropic() : Computes k-th derivative of first Barnett-Lothe integrand of an orthotropic medium
-  dkN2_orthotropic() : Computes k-th derivative of second Barnett-Lothe integrand of an R0-orthotropic medium
-  dkN1_r0_orthotropic() : Computes k-th derivative of first Barnett-Lothe integrand of an R0-orthotropic medium
-  dkN2_r0_orthotropic() : Computes k-th derivative of second Barnett-Lothe integrand of an R0-orthotropic medium
-  dkN1_square_symmetric() : Computes k-th derivative of first Barnett-Lothe integrand of an square-symmetric medium
-  dkN2_square_symmetric() : Computes k-th derivative of second Barnett-Lothe integrand of an R0-orthotropic medium
-  dkN1_isotropic() : Computes k-th derivative of first Barnett-Lothe integrand of an isotropic medium
-  dkN2_isotropic() : Computes k-th derivative of second Barnett-Lothe integrand of an isotropic medium
-  dh()	: Computes k-th derivative of angle dependent part of n-th gradient of Green's function
-  h() : Computes angle dependent part of k-th gradient of Green's function
-  dkh1i() : Computes k-th derivative of angle dependent part of 1-st gradient the Green's function
-  dknvec() : Computes k-th derivative of Cartesian component i of unit vector nvec
-  dkmvec() : Computes k-th derivative of Cartesian component i of unit vector mvec,  
   i.e. active counter-clockwise Pi/2 rotation of nvec
-  trapzd() : Performs numerical integration of Barnett-Lothe integrand by trapezoidal rule
-  qtrap() : Handles the trapezoidal rule for targeted error, i.e. 10.^-10

#### Prior version history:

 -  0.90 (May 31st, 2017),
 -  0.91 (June 7th, 2017).