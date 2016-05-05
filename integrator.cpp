// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "physics.hpp"
#include "integrator.hpp"

// -------------------------------------

using namespace Eigen;

// -------------------------------------


int F(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {

	flmatrix f = *(flmatrix *)fdata;

	Matrix7d I = Matrix7d::Zero(7,7);

	double k,tK,pK,pref;

	if (ndim == 2) {
		k = 1;
		tK = x[0];
		pK = x[1];
		pref = sin(tK)*fourierCoeff/(3-2*nKolmogorov);
	} else if (ndim == 3) {
		tK = x[1];
		pK = x[2];
		k = f.computeKfromKA(x[0]);
		pref = sin(tK)*fourierCoeff;
		if (k <= f.transK) {
			pref /= (3-2*nKolmogorov);
		} else {
			pref /= (3-2*nMHD);
		}
	}

	f.set_k(k,tK,pK);

	Matrix57d transform = MatrixXd::Zero(5,7);
	for (int j=0;j<3;j++) {
		transform(0,j) = f.a[j];
		transform(1,j) = f.b[j];
		transform(3,4+j) = f.a[j];
		transform(4,4+j) = f.b[j]; // Check if the velocity handling needs a c term as well
	}
	transform(2,3) = 1; // Density perturbation doesn't change under rotation

	I = pref*transform.transpose()*f.correlator*transform;

	int counter = 0;
	for (int j=0;j<7;j++) {
		for (int q=0;q<7;q++) {
			fval[counter] = I(j,q);
			counter++;
		}
	}

	return 0;
}