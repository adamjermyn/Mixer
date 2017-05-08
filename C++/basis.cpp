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

//	cout << db[0] << " " << db[1] << " " << db[2] << endl;

/*	
	cout << a[0] << " " << a[1] << " " << a[2] << endl;
	cout << b[0] << " " << b[1] << " " << b[2] << endl;
	cout << c[0] << " " << c[1] << " " << c[2] << endl;
	cout << dk[0] << " " << dk[1] << " " << dk[2] << endl;
	cout << db[0] << " " << db[1] << " " << db[2] << endl;
	cout << dc[0] << " " << dc[1] << " " << dd[2] << endl;
	cout << dd[0] << " " << dd[1] << " " << dc[2] << endl;
	cout << de[0] << " " << de[1] << " " << de[2] << endl;
*/

}
