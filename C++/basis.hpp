// -------------------------------------
// Header guard

#ifndef basis_h
#define basis_h

// -------------------------------------
// Includes

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"

// -------------------------------------

using namespace Eigen;

// -------------------------------------

/*

A basis object provides access to the basis
elements a-hat, b-hat, and c-hat as detailed
in our paper [insert citation]. These are provided
as normalised unit vectors.


This requires knowledge of the vector w-hat, which is specified
by the spherical polar angle, as well as k-hat, which
is specified by both the polar angle and the azimuthal
one.

The additional vectors d-hat and e-hat are defined as

d-hat = z-hat x c-hat
e-hat = z-hat x b-hat

These are not unit vectors!

Finally, the basis object provides a vector db, defined as

db = (1/|R Grad Omega|) d b-hat/dt,

and likewise

dc = (1/|R Grad Omega|) d c-hat/dt
dd = (1/|R Grad Omega|) d d-hat/dt
de = (1/|R Grad Omega|) d e-hat/dt

Note that a-hat is a constant in time so da is
not provided.

*/

// -------------------------------------

class basis
{

public:
	double kHat[3];

	double a[3];
	double b[3];
	double c[3];

	double d[3];
	double e[3];

	double db[3];
	double dc[3];
	double dd[3];
	double de[3];

	basis(tw);

	void set_k(double kT, double kP);

}