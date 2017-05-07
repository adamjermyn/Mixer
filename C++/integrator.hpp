// -------------------------------------
// Header guard

#include "cubature.h"

#ifndef integrator_h
#define integrator_h

// -------------------------------------
// Constants

const double pi2 = 6.28318530718;
const double fourierCoeff = 1;			 	// Accounts for the integrals being evaluated
							// in Fourier space. Careful derivation shows this to be
							// 1, but it is most easily seen by considering the
							// form of the delta function in a finite volume.
const double unitConv = 1/(pi2);			// This puts things in units of the mixing length.
							// The integration internally uses units of the mixing
							// wavevector, meaning that the lengths are all in units
							// of 1/k = L_mix/2pi. To fix this, we divide the position and
							// velocity correlators by 2pi to the power of the number of factors
							// of length, so they are in units of L_mix.

// -------------------------------------
// Functions

int integral(unsigned fdim, integrand f, void *fdata,
              unsigned dim, const double *xmin, const double *xmax, 
              size_t maxEval, double reqAbsError, double reqRelError, 
              error_norm norm,
              double *val, double *err);

int F(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

// -------------------------------------
// End Header guard
#endif
