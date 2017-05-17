// -------------------------------------
// Header guard

#ifndef physics_h
#define physics_h

// -------------------------------------
// Includes

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.hpp"
#include "basis.hpp"

// -------------------------------------

using namespace Eigen;

// -------------------------------------
// Constants

const double inf = 1.0/0.0;
const double pi = 3.141592653589;

// -------------------------------------
// Spectral indices (v ~ k^(-n))

const double nKolmogorov = 11./6;
const double powN = 3-2*nKolmogorov;

// -------------------------------------
// Adiabatic constant

const double gamma_ad = 5./3;

// -------------------------------------

class flmatrix
{

public:

	// Geometry
	basis ba = basis(0);

	double kmag;
	double kva;

	/*
	The angle theta defaults to pi/2 so that spherical and cylindrical
	coordinates are identical. Changing this to another value causes the
	returned coefficients to be in spherical coordinates evaluated at that
	polar angle.
	*/
	double theta = pi/2;

	// Physical parameters
	double va[3];
	double entHat[3];
	double presHat[3];
	double wmag;
	double N2;
	double omega;
	double transK;

	// Temporary arrays
	Matrix1 m = Matrix1::Zero();
	MatrixXcd eigvals;
	MatrixXcd eigvecs;
	MatrixXcd proj;

	// Intended outputs
	int* output; // 1 to output the given component, 0 otherwise.
	Matrix1 correlator = Matrix1::Zero();

	// Eigensolver
	EigenSolver<Matrix1> es;

	// Constructor
	flmatrix(double B, double tB, double pB, double w, double tW, double tS, double tP
		, double N22, double omegaa);

	void set_vecs();

	void set_M();

	Matrix1 derivative(int i);

	void compute_eigensystem();

	void compute_correlator();

	double computeKfromKA(double ka);

	// Set wavevector in spherical coordinates
	void set_k(double kmagg, double kT, double kP);

};

// -------------------------------------
// End Header guard
#endif
