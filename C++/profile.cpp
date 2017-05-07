// -------------------------------------

#include <iostream>
#include <fstream>

#include <stdlib.h>     /* strtod */

#define EIGEN_USE_BLAS

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "physics.hpp"

#include "integrator.hpp"

#include "cubature.h"

// -------------------------------------

using namespace Eigen;
using namespace std;

int main() {
	double tolr = 1e-11;
	double tola = 1e-11;
	double B = 0;
	double tB = 0;
	double pB = 0;
	double omega = 1.;
	double w = 1.;
	double tW = 0;
	double tS = 1.57;
	double tP = 1.57;
	double N2 = -1;
	int maxEval = 100000;

	const unsigned intDim = 2;
	const double transform_aa[3] = {0,0};
	const double transform_bb[3] = {pi,2*pi};

	double ret[36];
	double err[36];
	double* net = new double[72];
	
	flmatrix f(B,tB,pB,w,tW,tS,tP,N2,omega);
	integral(36,&F,&f,intDim,transform_aa,transform_bb,maxEval, tolr, tola ,ERROR_INDIVIDUAL,
			ret,err);

	for (int i=0;i<36;i++) {
		net[2*i] = ret[i];
		net[2*i + 1] = err[i];
		cout << ret[i] << " " << err[i] << endl;
	}

	return 0;
}

