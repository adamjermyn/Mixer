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

VectorC normalizeV(VectorC v, double eps) {
	//TODO: Add correction from c-hat piece
	double net = sqrt(eps+pow(abs(v(dim-1,0)),2) + pow(abs(v(dim-2,0)),2));	
	VectorC ret(v);
	v /= net;
	return v;
}

Matrix2 nullProjector(Matrix2 m, double eps) {
	/*
	This method takes as input a (dim*2 x dim*2) matrix and a threshold and returns
	the matrix which projects into its right null space,
	with the threshold used to distinguish null eigenvalues.
	*/

	JacobiSVD<Matrix2> svd(m, ComputeFullU | ComputeFullV );
	svd.compute(m);

	Vector2 vals = Vector2::Zero();
	Matrix2 ret = Matrix2::Zero();

	for (int i=0;i<10;i++) {
		if (abs(svd.singularValues()(i)) < eps) {
			ret += ((svd.matrixV().col(i))*(svd.matrixV().col(i)).adjoint()).real();
		}
	}

	return ret;
}

double vGrowth(VectorC2 v) {
	/*
	This method takes as input a vector of size 10 and returns the growth rate
	implied by the elements of the vector, assuming that the first five are
	(alpha, beta, d log rho, alpha-dot, beta-dot) and the next five are the
	time derivatives of these.
	*/
	cdouble g0 = v[3]*conj(v[8]) + v[4]*conj(v[9]);
	cdouble g1 = conj(g0);
	return ((g0 + g1)/(eps + abs(v[3])*abs(v[3]) + abs(v[4])*abs(v[4]))).real();

}