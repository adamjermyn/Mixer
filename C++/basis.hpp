// -------------------------------------
// Header guard

#ifndef basis_h
#define basis_h

// -------------------------------------
// Includes

#include "linalg.hpp"

// -------------------------------------

/*

A basis object provides access to the basis
elements a-hat, b-hat, and c-hat as detailed
in Jermyn et al 2017. These are provided
as normalised unit vectors.


This requires knowledge of the vector w-hat, which is specified
by the spherical polar angle, as well as k-hat, which
is specified by both the polar angle and the azimuthal
one.

The additional vectors d-hat and e-hat are defined as

e-hat = z-hat x b-hat
d-hat = z-hat x db

These are not unit vectors!

The basis object provides a vector dk, defined as

dk = (1/|R Grad Omega|) d k-hat/dt

and likewise

db = (1/|R Grad Omega|) d b-hat/dt,
dd = (1/|R Grad Omega|) d d-hat/dt
de = (1/|R Grad Omega|) d e-hat/dt

Note that a-hat is a constant in time so da is
not provided.

Finally, we store the inner product of kHat and wHat as kw
and the inner product of dk and wHat as dkw.

*/

// -------------------------------------

const int maxOrder = 1;
const double zhat[3] = {0,0,1};

// -------------------------------------

class basis
{

public:
	double wHat[3];
	double kHat[3];

	double a[3];
	double b[3];
	double c[3];

	double d[3];
	double e[3];

	double dk[maxOrder + 1][3];
	double db[maxOrder + 1][3];
	double de[maxOrder + 1][3];

	double dd[maxOrder][3];

	double tW;

	basis(double tWW);

	void set_k(double kT, double kP);

};

// -------------------------------------
// End Header guard
#endif
