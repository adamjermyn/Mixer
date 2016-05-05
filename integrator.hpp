// -------------------------------------
// Header guard

#ifndef integrator_h
#define integrator_h

// -------------------------------------
// Constants

const double pi = 3.141592653589;
const double fourierCoeff = 1/pow(2*pi,3);

// -------------------------------------
// Functions

int F(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

// -------------------------------------
// End Header guard
#endif