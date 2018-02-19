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

int main() {
	double tolr = 1e-3;
	double tola = 1e-12;
	double B = 0;
	double tB = 0;
	double pB = 0;
	double omega = 0.1;
	double w = 0.01;
	double tW = 0.75;
	double tS = 1.57;
	double tP = 1.57;
	double N2 = -1;
	int maxEval = 1000000;

	const unsigned intDim = 2;
	const double transform_aa[2] = {0,0};
	const double transform_bb[2] = {pi,2*pi};

	double ret[36];
	double err[36];

	int output[36];
	for (int i=0;i<36;i++) {
		output[i] = 1;
	}
	
	flmatrix f(B,tB,pB,w,tW,tS,tP,N2,omega,1);
	f.theta = 1.57;
	f.output = output;
	integral(36,&F,&f,intDim,transform_aa,transform_bb,maxEval, tolr, tola ,ERROR_INDIVIDUAL,
			ret,err);

//	for (int i=0;i<36;i++) {
//		cout << ret[i] << " " << err[i] << endl;
//	}

	return 0;
}

