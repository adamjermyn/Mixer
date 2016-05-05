// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "physics.hpp"

#include <iostream>
#include "TasmanianSparseGrid.hpp"


// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

Matrix7d integrate(flmatrix f, int intDim, int num_points, double* points, double* weights) {
	Matrix7d I = Matrix7d::Zero(7,7);
	double k=1;
	for (int i=0;i<num_points;i++) {
		double tK = points[i*intDim];
		double pK = points[i*intDim+1];
//		double k = points[i*intDim+2];

		f.set_k(k,tK,pK); // TODO: Make the integration go over the proper range (k=1 to infinity)

		Matrix57d transform = MatrixXd::Zero(5,7);
		for (int j=0;j<3;j++) {
			transform(0,j) = f.a[j];
			transform(1,j) = f.b[j];
			transform(3,4+j) = f.a[j];
			transform(4,4+j) = f.b[j]; // Check if the velocity handling needs a c term as well
		}
		transform(2,3) = 1; // Density perturbation doesn't change under rotation

		I += weights[i]*sin(tK)*transform.transpose()*f.correlator*transform;

	}
	return I;
}

int F(unsigned ndim, const double *x, void *fdata,
      unsigned fdim, double *fval) {

	flmatrix f = *(flmatrix *)fdata;
//	static_cast<flmatrix>(*fdata);

	Matrix7d I = Matrix7d::Zero(7,7);
	double k=1;
	double tK = x[0];
	double pK = x[1];
//		double k = points[i*intDim+2];

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