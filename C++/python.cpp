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

extern "C" double* coeffs3(double tolr, double tola, double B, double tB, double pB, double omega, double w, double tW, double tS, double tP, double N2, int maxEval) {

	const unsigned intDim = 3;
	const double transform_aa[3] = {0,0,0};
	const double transform_bb[3] = {1,pi,2*pi};

	double ret[correlDim*correlDim];
	double err[correlDim*correlDim];
	double* net = new double[2*correlDim*correlDim];
	
	flmatrix f(B,tB,pB,w,tW,tS,tP,N2,omega);
	integral(correlDim*correlDim,&F,&f,intDim,transform_aa,transform_bb,maxEval, tolr, tola ,ERROR_INDIVIDUAL,
			ret,err);

	for (int i=0;i<correlDim*correlDim;i++) {
		net[2*i] = ret[i];
		net[2*i + 1] = err[i];
	}

	return net;
}

extern "C" double* coeffs3spherical(double theta, double tolr, double tola, double B, double tB, double pB, double omega, double w, double tW, double tS, double tP, double N2, int maxEval) {

	const unsigned intDim = 3;
	const double transform_aa[3] = {0,0,0};
	const double transform_bb[3] = {1,pi,2*pi};

	double ret[correlDim*correlDim];
	double err[correlDim*correlDim];
	double* net = new double[2*correlDim*correlDim];
	
	flmatrix f(B,tB,pB,w,tW,tS,tP,N2,omega);
	f.theta = theta;
	integral(correlDim*correlDim,&F,&f,intDim,transform_aa,transform_bb,maxEval, tolr, tola ,ERROR_INDIVIDUAL,
			ret,err);

	for (int i=0;i<correlDim*correlDim;i++) {
		net[2*i] = ret[i];
		net[2*i + 1] = err[i];
	}

	return net;
}

extern "C" double* coeffs2(double tolr, double tola, double omega, double w, double tW, double tS, double tP, double N2, int maxEval) {

	double B = 0;
	double tB = 0;
	double pB = 0;

	const unsigned intDim = 2;
	const double transform_aa[3] = {0,0};
	const double transform_bb[3] = {pi,2*pi};

	double ret[correlDim*correlDim];
	double err[correlDim*correlDim];
	double* net = new double[2*correlDim*correlDim];
	
	flmatrix f(B,tB,pB,w,tW,tS,tP,N2,omega);
	integral(correlDim*correlDim,&F,&f,intDim,transform_aa,transform_bb,maxEval, tolr, tola,ERROR_INDIVIDUAL,
			ret,err);

	for (int i=0;i<correlDim*correlDim;i++) {
		net[2*i] = ret[i];
		net[2*i + 1] = err[i];
	}

	return net;
}

extern "C" double* coeffs2spherical(double theta, double tolr, double tola, double omega, double w, double tW, double tS, double tP, double N2, int maxEval) {

	double B = 0;
	double tB = 0;
	double pB = 0;

	const unsigned intDim = 2;
	const double transform_aa[3] = {0,0};
	const double transform_bb[3] = {pi,2*pi};

	double ret[correlDim*correlDim];
	double err[correlDim*correlDim];
	double* net = new double[2*correlDim*correlDim];
	
	flmatrix f(B,tB,pB,w,tW,tS,tP,N2,omega);
	f.theta = theta
	integral(correlDim*correlDim,&F,&f,intDim,transform_aa,transform_bb,maxEval, tolr, tola,ERROR_INDIVIDUAL,
			ret,err);

	for (int i=0;i<correlDim*correlDim;i++) {
		net[2*i] = ret[i];
		net[2*i + 1] = err[i];
	}

	return net;
}
