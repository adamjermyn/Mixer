// -------------------------------------

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "basis.hpp"

#include "linalg.hpp"

#include <iostream>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_fact.h>

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

	cross(kTemp, a, db[0]);

	for (int i=1;i<=maxOrder+2;i++) {
		dk[i][0] = gsl_sf_fact(i) * cos(p) * gsl_sf_legendre_Pl(i, -cos(t)) * sin(t);
		dk[i][1] = gsl_sf_fact(i) * sin(p) * gsl_sf_legendre_Pl(i, -cos(t)) * sin(t);
		dk[i][2] = gsl_sf_fact(i) * (cos(t) * gsl_sf_legendre_Pl(i, -cos(t)) + gsl_sf_legendrePl(i-1, -cos(t)));

		db[i][0] = cos(p) * (gsl_sf_legendre_Pl(i-1, -cos(t)) + cos(t) * gsl_sf_legendre_Pl(i, -cos(t)));
		db[i][1] = sin(p) * (gsl_sf_legendre_Pl(i-1, -cos(t)) + cos(t) * gsl_sf_legendre_Pl(i, -cos(t)));
		db[i][2] = -sin(t) * gsl_sf_legendre_Pl(i, -cos(t));
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
	cross(zhat,b,e);

	// Now we can just use the definitions of these vectors.
	for(int i=0;i<maxOrder + 2;i++) {
		cross(zhat, db[i], de[i]);
	}

}
