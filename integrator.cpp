// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "physics.hpp"

#include <iostream>


// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

int F(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {

	flmatrix f = *(flmatrix *)fdata;

	Matrix7d I = Matrix7d::Zero(7,7);

	double k,tK,tP;

	if (ndim == 2) {
		k = 1;
		tK = x[0];
		tP = x[1];
	} else if (ndim == 3) {
		k = x[0];
		tK = x[1];
		tP = x[2];
	}

	f.set_k(k,tK,pK); // TODO: Make the integration go over the proper range (k=1 to infinity)

	Matrix57d transform = MatrixXd::Zero(5,7);
	for (int j=0;j<3;j++) {
		transform(0,j) = f.a[j];
		transform(1,j) = f.b[j];
		transform(3,4+j) = f.a[j];
		transform(4,4+j) = f.b[j]; // Check if the velocity handling needs a c term as well
	}
	transform(2,3) = 1; // Density perturbation doesn't change under rotation

	I = sin(tK)*transform.transpose()*f.correlator*transform;

	int counter = 0;
	for (int j=0;j<7;j++) {
		for (int q=0;q<7;q++) {
			fval[counter] = I(j,q);
			counter++;
		}
	}

	return 0;
}