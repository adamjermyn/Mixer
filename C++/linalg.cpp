// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "linalg.hpp"

// -------------------------------------

using namespace Eigen;

// -------------------------------------

double dot(const double a[3], const double b[3]) {
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

void cross(const double a[3], const double b[3], double ret[3]) {
	ret[0] = a[1]*b[2]-a[2]*b[1];
	ret[1] = a[2]*b[0]-a[0]*b[2];
	ret[2] = a[0]*b[1]-a[1]*b[0];
}

void normalize(double* v, double eps) {
	double norm = sqrt(dot(v,v)+eps);
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;
}

void sphericalToCartesian(double r, double t, double p, double ret[3]) {
	ret[0] = r*sin(t)*cos(p);
	ret[1] = r*sin(t)*sin(p);
	ret[2] = r*cos(t);
}

Vector5cd normalizeV(Vector5cd v, double eps) {
	double net = sqrt(eps+pow(abs(v(dim-1,0)),2) + pow(abs(v(dim-2,0)),2));	
	Vector5cd ret(v);
	v /= net;
	return v;
}

Matrix5cd normalizeM(Matrix5cd m, double eps) {
	// This function takes as input a matrix with eigenvectors as columns
	// and returns a copy with each column normalized such that
	// the sum of the norm squares of the last two elements is unity.
	Matrix5cd ret(m);
	double net;
	double elem;
	//We've optimized dim for matrices of size five
	for (int i=0;i<dim;i++) {
		net = sqrt(eps+pow(abs(m(dim-1,i)),2) + pow(abs(m(dim-2,i)),2));
		ret.col(i) /= net;
	}
	return ret;
}


Matrix10d nullProjector(Matrix10d m, double eps) {
	/*
	This method takes as input a 10x10 matrix and a threshold and returns
	the matrix which projects into the null space of the given matrix,
	with the threshold used to distinguish null eigenvalues.
	*/

	EigenSolver<Matrix10d> es10;
	es.compute(m);

	Matrix10d ret = Matrix10d::Zero();

	for (int i=0;i<10;i++) {
		if (abs(es.eigenvalues()(i)) < eps) {
			ret += (es.eigenvectors().col(i))*(es.eigenvectors().col(i)).adjoint();
		}
	}

	return ret;
}
