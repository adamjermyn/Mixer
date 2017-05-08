// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "linalg.hpp"

#include <iostream>

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

double dot(const double a[3], const double b[3]) {
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

void cross(const double a[3], const double b[3], double ret[3]) {
	ret[0] = a[1]*b[2]-a[2]*b[1];
	ret[1] = a[2]*b[0]-a[0]*b[2];
	ret[2] = a[0]*b[1]-a[1]*b[0];
}

void normalize(double* v) {
	double norm = vectorNormEPS + sqrt(dot(v,v));
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;
}

void sphericalToCartesian(double r, double t, double p, double ret[3]) {
	ret[0] = r*sin(t)*cos(p);
	ret[1] = r*sin(t)*sin(p);
	ret[2] = r*cos(t);
}

VectorC normalizeV(VectorC v, double kPhi, double w) {
	/*
	This method returns a normalized version of the given state vector,
	such that the velocity is unity. The velocity is taken to reside
	in the final two elements of the vector.

	The normalization is computed as

	|v|^2 = |v[dim - 2]|^2 + |v[dim - 1]|^2  + |v[1]|^2 * k_phi^2 w^2

	*/

	double net = pow(abs(v(dim - 2)), 2) + pow(abs(v(dim - 1)), 2) + pow(abs(v(1)), 2)*pow(kPhi*w, 2);
	net = sqrt(net);
	VectorC ret(v);
	ret /= (velocityNormEPS + net);
	return ret;
}

Eigen::MatrixXcd nullProjector(Matrix2 m) {
	/*
	This method takes as input a (dim*2 x dim*2) matrix and a threshold and returns
	the matrix which projects into its right null space,
	with the threshold used to distinguish null eigenvalues.
	*/

	JacobiSVD<Matrix2> svd(m, ComputeFullU | ComputeFullV );
	svd.compute(m);

	// Count null singular values
	int num = 0;
	for (int i=0;i<2*dim;i++) {
		if (abs(svd.singularValues()(i)) < svdEPS) {
			num += 1;
		}
	}

	MatrixXcd temp = MatrixXcd::Zero(num, 2*dim);

	// Construct projector
	int counter = 0;
	for (int i=0;i<2*dim;i++) {
		if (abs(svd.singularValues()(i)) < svdEPS) {
			for (int j=0;j<2*dim;j++) {
				temp(counter,j) += svd.matrixV().adjoint()(i,j);
			}
			counter += 1;
		}
	}

	return temp;

}
