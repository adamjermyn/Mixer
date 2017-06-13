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

void rotY(double theta, double v[3]) {
	// Rotates the specified vector an angle theta about the middle axis
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

VectorC normalizeV(VectorC v, double kPhi, double w, double kw, double eps) {
	/*
	This method returns a normalized version of the given state vector,
	such that the velocity is unity. The velocity is taken to reside
	in the final two elements of the vector.

	The normalization is computed as

	|v|^2 = |v[dim - 2]|^2 + |v[dim - 1]|^2  + |v[1]|^2 * k_phi^2 w^2 + 2*Real[v[1]*v[dim-1]*kPhi*w*kw]

	where the final term arises because b-hat is not perpendicular to c-hat.

	*/

	double net = pow(abs(v(dim - 2)), 2) + pow(abs(v(dim - 1)), 2) + pow(abs(v(1)), 2)*pow(kPhi*w, 2) + 2*(v(dim - 1)*v(1)*kPhi*w*kw).real();
	net = sqrt(eps + net);
	VectorC ret(v);
	ret /= net;
	return ret;
}

