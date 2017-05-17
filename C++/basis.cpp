// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "basis.hpp"

#include "linalg.hpp"

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
	double t = acos(kTemp[2]);
	double p = atan2(kTemp[1], kTemp[0]);

	// Now we apply the analytic normalized cross product
	a[0] = sin(p);
	a[1] = -cos(p);
	a[2] = 0;

	dk[0][0] = kTemp[0];
	dk[0][1] = kTemp[1];
	dk[0][2] = kTemp[2];
	if (maxOrder > 0) {
		dk[1][0] = -cos(p)*cos(t)*sin(t);
		dk[1][1] = -cos(t)*sin(p)*sin(t);
		dk[1][2] = sin(t)*sin(t);
	}
	if (maxOrder > 1) {
		dk[2][0] = 0.5*cos(p)*(1+3*cos(2*t))*sin(t);
		dk[2][1] = 0.5*(1+3*cos(2*t))*sin(p)*sin(t);
		dk[2][2] = -3*cos(t)*sin(t)*sin(t);
	}
	if (maxOrder > 2) {
		dk[3][0] = -1.5*cos(p)*cos(t)*(5*cos(2*t) - 1)*sin(t);
		dk[3][1] = -1.5*cos(t)*(5*cos(2*t) - 1)*sin(p)*sin(t);
		dk[3][2] = 1.5*(5*cos(2*t) + 3)*sin(t)*sin(t);
	}

	cross(kHat, a, db[0]);


	if (maxOrder > 0) {
		db[1][0] = cos(p)*pow(sin(t),2);
		db[1][1] = sin(p)*pow(sin(t),2);
		db[1][2] = cos(t)*sin(t);
	}
	if (maxOrder > 1) {
		db[2][0] = -3*cos(p)*cos(t)*pow(sin(t),2);
		db[2][1] = -3*sin(p)*cos(t)*pow(sin(t),2);
		db[2][2] = (1-3*pow(cos(t),2))*sin(t);
	}
	if (maxOrder > 2) {	
		db[3][0] = 1.5*cos(p)*(3+5*cos(2*t))*pow(sin(t),2);
		db[3][1] = 1.5*sin(p)*(3+5*cos(2*t))*pow(sin(t),2);
		db[3][2] = 3*cos(t)*sin(t)*(-3+5*pow(cos(t),2));
	}

	// And finally we counter-rotate
	rotY(-tW, a);
	for (int i=0;i<maxOrder + 1;i++) {
		rotY(-tW, dk[i]);
		rotY(-tW, db[i]);
	}

	// Now a is orthogonal to both kHat and wHat, so
	// the results of subsequent cross-products don't
	// require normalization.
	cross(kHat,a,b);
	cross(wHat,a,c);

	// Helper vectors
	cross(zhat,c,d);
	cross(zhat,b,e);

	// Now we can just use the definitions of these vectors.
	// Note that it matters that the normalization factor on db
	// equals that on dc (i.e. that wHat is orthogonal to b).
	for(int i=0;i<maxOrder + 1;i++) {
		cross(wHat, db[i], dc[i]);
		cross(zhat, dc[i], dd[i]);
		cross(zhat, db[i], de[i]);
	}

}
