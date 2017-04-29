// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "basis.hpp"

#include "linalg.hpp"

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

basis(double tw) {
	sphericalToCartesian(1,tW,0,wHat);
}

void set_k(double kT, double kP) {

	// Construct k-hat
	sphericalToCartesian(1,kT,kP,kHat);

	// Construct basis
	cross(kHat,wHat,a);
	cross(kHat,a,b);
	cross(wHat,a,c);

	// Normalize basis
	normalize(a,eps);
	normalize(b,eps);
	normalize(c,eps);	

	// The original version had these cross products
	// coming before the normalization, so check this!

	// Helper vectors
	cross(zhat,c,d);
	cross(zhat,b,e);

	// Derivative vectors
	// Finish writing this
}