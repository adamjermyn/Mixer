// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "basis.hpp"

#include "linalg.hpp"

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

basis::basis(double tW) {
	sphericalToCartesian(1,tW,0,wHat);
}

void basis::set_k(double kT, double kP) {

	// Construct k-hat
	sphericalToCartesian(1,kT,kP,kHat);

	// Construct basis
	cross(kHat,wHat,a);
	cross(kHat,a,b);
	cross(wHat,a,c);

	// Normalize basis
	normalize(a);
	normalize(b);
	normalize(c);	

	// The original version had these cross products
	// coming before the normalization, so check this!

	// Helper vectors
	cross(zhat,c,d);
	cross(zhat,b,e);

	// Derivative vectors
	kw = dot(kHat, wHat);

	// dk
	for (int i=0;i<3;i++)
		dk[i] = kHat[1] * (wHat[i] + kHat[i]*kw);

	dkw = dot(dk, wHat);
	double denom = eps + pow(1 - kw*kw,0.5);

	// db
	for (int i=0;i<3;i++) {
		db[i] = dk[i]*kw/denom;
		db[i] += kHat[i]*dkw/denom;
		db[i] += b[i]*dkw*kw/pow(denom, 2);
	}

	// Now we can just use the definitions of these vectors.
	// Note that it matters that the normalization factor on db
	// equals that on dc (i.e. that wHat is orthogonal to b).
	cross(wHat, db, dc);
	cross(zhat, dc, dd);
	cross(zhat, db, de);

}
