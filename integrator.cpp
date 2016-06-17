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

	double maxErrH = 0;

	hcubature(fdim,f,fdata,dim,xmin,xmax,maxEval,reqAbsError,reqRelError,norm,val,err);

	for (int i=0;i<fdim;i++) {
		if (err[i] > maxErrH)
			maxErrH = err[i];
	}

	double maxErrP = 0;

	pcubature(fdim,f,fdata,dim,xmin,xmax,maxEval,reqAbsError,reqRelError,norm,val,err);

	for (int i=0;i<fdim;i++) {
		if (err[i] > maxErrP)
			maxErrP = err[i];
	}

	if (maxErrP < maxErrH) {
		pcubature(fdim,f,fdata,dim,xmin,xmax,0,reqAbsError,reqRelError,norm,val,err);
	} else {
		hcubature(fdim,f,fdata,dim,xmin,xmax,0,reqAbsError,reqRelError,norm,val,err);
	}

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

	// We now bound the position correlation functions to be at most 1.
	// This means that the position correlations are never over more than
	// a single mixing length (in practice correlations which exceed this length
	// will be distorted by longer-scale processes, making our calculations invalid
	// if we don't do this). The max is because failure of perturbation theory
	// can cause negative terms on the diagonal which are unphysical, so we don't
	// want these to impact the correction. These terms arise when the matrix is nearly
	// degenerate, in which case the growth correction makes them small relative to the
	// other contributions to the diagonal upon integration.

	double net = sqrt(1 + fmax(0,I(0,0) + I(1,1) + I(2,2)));

	I.row(0) /= net;
	I.row(1) /= net;
	I.row(2) /= net;
	I.col(0) /= net;
	I.col(1) /= net;
	I.col(2) /= net;

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