// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "linalg.hpp"

#include <iostream>

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

int nCr(int n, int k) {
	if ((k > n) || (k < 0)) {
		return 0;
	} else if (k == n) {
		return 1;
	} else {
		return factorial(n)/(factorial(n-k)*factorial(k));		
	}
}

double dot(const double a[3], const double b[3]) {
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

cdouble dot(const cdouble a[3], const cdouble b[3]) {
	return conj(a[0])*b[0]+conj(a[1])*b[1]+conj(a[2])*b[2];
}

void rotY(double theta, double v[3]) {
	// Rotates the specified vector an angle minus theta about the middle axis
	double x, y;
	x = v[0]*cos(theta) + v[2]*sin(theta);
	y = -v[0]*sin(theta) + v[2]*cos(theta);
	v[0] = x;
	v[2] = y;
}

void cross(const double a[3], const double b[3], double ret[3]) {
	ret[0] = a[1]*b[2]-a[2]*b[1];
	ret[1] = a[2]*b[0]-a[0]*b[2];
	ret[2] = a[0]*b[1]-a[1]*b[0];
}

void sphericalToCartesian(double r, double t, double p, double ret[3]) {
	ret[0] = r*sin(t)*cos(p);
	ret[1] = r*sin(t)*sin(p);
	ret[2] = r*cos(t);
}

void conjugate(cdouble v[3]) {
	v[0] = conj(v[0]);
	v[1] = conj(v[1]);
	v[2] = conj(v[2]);
}
