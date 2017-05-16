// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "basis.hpp"

#include "linalg.hpp"
#include "mathematica.hpp"

#include <iostream>

// -------------------------------------

using namespace Eigen;
using namespace std;

// -------------------------------------

basis::basis(double tWW) {
	tW = tWW;
	sphericalToCartesian(1,tW,0,wHat);
}

void basis::set_k(double kT, double kP) {

	// Construct k-hat
	sphericalToCartesian(1,kT,kP,kHat);

	// Construct a = normalize(cross(k,w))
	
	// First we construct a temporary copy of kHat
	double kTemp[3];
	kTemp[0] = kHat[0];
	kTemp[1] = kHat[1];
	kTemp[2] = kHat[2];

	// First we rotate the coordinate system about the middle
	// axis so that wHat = zHat.
	rotY(tW, kTemp);

	// Next we extract the spherical angles of kTemp in this new
	// coordinate system.
	double p = atan2(kTemp[0], kTemp[2]);
	double t = atan2(sqrt(kTemp[0]*kTemp[0] + kTemp[1]*kTemp[1]), kTemp[2]);

	// Now we apply the analytic normalized cross product
	a[0] = sin(p);
	a[1] = -cos(p);
	a[2] = 0;

	// And finally we counter-rotate
	rotY(-tW, a);

	// Now a is orthogonal to both kHat and wHat, so
	// the results of subsequent cross-products don't
	// require normalization.
	cross(kHat,a,b);
	cross(wHat,a,c);

	// Helper vectors
	cross(zhat,c,d);
	cross(zhat,b,e);

	// Derivative vectors
	kw = dot(kHat, wHat);

	// dk
	for (int i=0;i<3;i++)
		dk[i] = kHat[1] * (wHat[i] + kHat[i]*kw);

	dkw = dot(dk, wHat);

	// db
	db[0] = db0(kT, kP, tW);
	db[1] = db1(kT, kP, tW);
	db[2] = db2(kT, kP, tW);

	// Now we can just use the definitions of these vectors.
	// Note that it matters that the normalization factor on db
	// equals that on dc (i.e. that wHat is orthogonal to b).
	cross(wHat, db, dc);
	cross(zhat, dc, dd);
	cross(zhat, db, de);

}
