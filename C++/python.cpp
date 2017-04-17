// -------------------------------------

#include <iostream>
#include <fstream>

#include <stdlib.h>     /* strtod */

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "physics.hpp"

#include "integrator.hpp"

#include "cubature.h"

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

const int maxEval = 100000;

// -------------------------------------

extern "C" double* coeffs3(double tolr, double tola, double B, double tB, double pB, double omega, double w, double tW, double tS, double tP, double N2, double chi) {

	const unsigned intDim = 3;
	const double transform_aa[3] = {0,0,0};
	const double transform_bb[3] = {1,pi,2*pi};

	double ret[49];
	double err[49];
	double* net = new double[98];
	
	flmatrix f(B,tB,pB,w,tW,tS,tP,N2,chi,omega);
	integral(49,&F,&f,intDim,transform_aa,transform_bb,maxEval,3e-4,3e-4,ERROR_INDIVIDUAL,
			ret,err);

	for (int i=0;i<49;i++) {
		net[2*i] = ret[i];
		net[2*i + 1] = err[i];
	}

	return net;
}

extern "C" double* coeffs2(double tolr, double tola, double omega, double w, double tW, double tS, double tP, double N2) {

	double B = 0;
	double tB = 0;
	double pB = 0;

	double chi = 0;

	const unsigned intDim = 2;
	const double transform_aa[3] = {0,0};
	const double transform_bb[3] = {pi,2*pi};

	double ret[49];
	double err[49];
	double* net = new double[98];
	
	flmatrix f(B,tB,pB,w,tW,tS,tP,N2,chi,omega);
	integral(49,&F,&f,intDim,transform_aa,transform_bb,maxEval,3e-4,3e-4,ERROR_INDIVIDUAL,
			ret,err);

	for (int i=0;i<49;i++) {
		net[2*i] = ret[i];
		net[2*i + 1] = err[i];
	}

	return net;
}
