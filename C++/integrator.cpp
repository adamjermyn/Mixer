// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "physics.hpp"
#include "integrator.hpp"

#include <iostream>

#include "cubature.h"

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

int integral(unsigned fdim, integrand f, void *fdata,
              unsigned dim, const double *xmin, const double *xmax, 
              size_t maxEval, double reqAbsError, double reqRelError, 
              error_norm norm,
              double *val, double *err) {

	hcubature(fdim,f,fdata,dim,xmin,xmax,maxEval,reqRelError,reqAbsError,norm,val,err);

	return 0;
}


int F(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {

	flmatrix f = *(flmatrix *)fdata;

	Matrix7d I = Matrix7d::Zero(7,7);

	double k,tK,pK,pref;

	if (ndim == 2) {
		k = 1;
		tK = x[0];
		pK = x[1];
		pref = -sin(tK)*fourierCoeff/(3-2*nKolmogorov);
	} else if (ndim == 3) {
		tK = x[1];
		pK = x[2];
		k = f.computeKfromKA(x[0]);
		pref = sin(tK)*fourierCoeff;
		if (k <= f.transK) {
			pref /= -(3-2*nKolmogorov);
		} else {
			pref /= -(3-2*nMHD);
		}
	}

	f.set_k(k,tK,pK);

	Matrix57d transform = MatrixXd::Zero(5,7);
	for (int j=0;j<3;j++) {
		transform(0,j) = f.a[j];
		transform(1,j) = f.b[j];
		transform(3,4+j) = f.a[j];
		transform(4,4+j) = f.b[j];
		transform(1,4+j) = f.c[j]*f.wmag*f.omega*f.kHat[1]; // The velocity has an additional term due to the sheared coordinate system.
	}
	transform(2,3) = 1; // Density perturbation doesn't change under rotation
	Matrix75d transformT = transform.transpose();

	// Apply transform to real-space coordinates

	I = transformT*f.correlator*transform;

	// Apply unit conversion

	I *= unitConv*unitConv; // Needs to have two factors because we have two length factors

	// Apply integration prefactor

	I *= pref;

	// Place output in the output array

	for (int j=0;j<7;j++) {
		for (int q=0;q<7;q++) {
			fval[7*j+q] = I(j,q);
		}
	}

	// The following corrects the units on the density correlator:

	for (int j=0;j<7;j++) {
		fval[7*j+3] /= unitConv;
		fval[7*3+j] /= unitConv;
	}

	return 0;
}
