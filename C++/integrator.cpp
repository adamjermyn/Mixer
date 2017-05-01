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

	MatrixCorr I = MatrixCorr::Zero();

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

	// Compute physics
	f.set_k(k,tK,pK);

	// Construct coordinate transform
	MatrixRegCorr transform = MatrixRegCorr::Zero();
	for (int j=0;j<3;j++) {
		transform(0,j) = f.ba.a[j];
		transform(1,j) = f.ba.b[j];
		transform(2,4+j) = f.ba.a[j];
		transform(3,4+j) = f.ba.b[j];
		transform(1,4+j) = f.ba.c[j]*f.wmag*f.omega*f.ba.kHat[1]; // The velocity has an additional term due to the sheared coordinate system.
	}
	transform(2,3) = 1; // Density perturbation doesn't change under rotation
	MatrixCorrReg transformT = transform.transpose();

	// Apply transform to real-space coordinates

	I = transformT*f.correlator*transform;

	// Apply unit conversion

	I *= unitConv*unitConv; // Needs to have two factors because we have two length factors

	// Apply integration prefactor

	I *= pref;

	// Place output in the output array

	for (int j=0;j<correlDim;j++) {
		for (int q=0;q<correlDim;q++) {
			fval[correlDim*j+q] = I(j,q);
		}
	}

	return 0;
}
