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
		pref = -sin(tK)*fourierCoeff/powN1;
	} else if (ndim == 3) {
		tK = x[1];
		pK = x[2];
		k = f.computeKfromKA(x[0]);
		pref = -sin(tK)*fourierCoeff;
		
		if (k <= f.transK) { 
			pref /= powN1;
		} else { 
			pref /= powN2;
		} 

	}

	// Compute physics
	f.set_k(k,tK,pK);

	// Construct coordinate transform
	MatrixRegCorr transform = MatrixRegCorr::Zero();
	for (int j=0;j<3;j++) {
		transform(0,j) = f.ba.a[j];
		transform(1,j) = f.ba.b[j];
		transform(2,3+j) = f.ba.a[j];
		transform(3,3+j) = f.ba.b[j];
		transform(1,3+j) = f.ba.c[j]*f.wmag*f.ba.kHat[1]; // The velocity has this additional term due to the sheared coordinate system.
	}

	MatrixCorrReg transformT = transform.transpose();

	// Apply transform to real-space coordinates

	I = transformT*f.correlator*transform;

	// Construct cylindrical to spherical transform

	MatrixSpace temp = MatrixSpace::Zero();
	MatrixSpace2 transformR = MatrixSpace2::Zero();

	double theta = f.theta;
	temp(0,0) = sin(theta);
	temp(0,2) = cos(theta);
	temp(1,0) = cos(theta);
	temp(1,2) = -sin(theta);
	temp(2,1) = 1;
	for (int i=0;i<spaceDim;i++) {
		for (int j=0;j<spaceDim;j++) {
			transformR(i,j) = temp(i,j);
			transformR(i + spaceDim, j + spaceDim) = temp(i,j);
		}
	}

	MatrixSpace2 transformTR = transformR.transpose();

	// Apply transform

	I = transformR * I * transformTR;

	// Apply unit conversion (needs to have two factors because we have two length factors)

	I *= unitConv*unitConv;

	// Apply integration prefactor

	I *= pref;

	// Place output in the output array

	for (int j=0;j<correlDim;j++) {
		for (int q=0;q<correlDim;q++) {
			if (f.output[correlDim*j + q] == 1) {
				fval[correlDim*j+q] = I(j,q);
			} else {
				fval[correlDim*j+q] = 0;
			}
		}
	}

	return 0;
}
