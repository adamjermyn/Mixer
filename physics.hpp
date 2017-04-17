// -------------------------------------
// Header guard

#ifndef physics_h
#define physics_h

// -------------------------------------
// Includes
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "linalg.hpp"

// -------------------------------------

using namespace Eigen;

// -------------------------------------
// Constants

const double zhat[3] = {0,0,1};
const double inf = 1.0/0.0;

// -------------------------------------
// Spectral indices (v ~ k^(-n))

const double nKolmogorov = 11./6;
const double nMHD = 8./3;
const double pow1 = 3-2*nKolmogorov;
const double pow2 = 3-2*nMHD;

// -------------------------------------
// Adiabatic constant

const double gamma = 5./3;

// -------------------------------------

class flmatrix
{

public:

		// Geometry
	double k[3];

	double kmag;
	double kHat[3];

	double a[3];
	double b[3];
	double c[3];
	double d[3];
	double e[3];

	double kva;


	// Physical parameters
	double va[3];
	double wHat[3];
	double entHat[3];
	double presHat[3];
	double wmag;
	double N2;
	double chi;
	double omega;
	double transK;

	// Intended outputs
	Matrix5d m;
	Matrix5d mdot;
	Matrix10cd eigvals;
	Matrix10cd eigvecs;
	Matrix5cd proj;
	Matrix5d correlator;

	// Eigensolver
	EigenSolver<Matrix5d> es;
	EigenSolver<Matrix10d> es10;

	// Constructor
	flmatrix(double B, double tB, double pB, double w, double tW, double tS, double tP
		, double N22, double chii, double omegaa);

	void set_vecs();

	void set_M();

	void set_Mdot();

	void compute_eigensystem();

	void compute_correlator();

	double computeKfromKA(double ka);

	// Set wavevector in spherical coordinates
	void set_k(double kmagg, double kT, double kP);

};

// -------------------------------------
// End Header guard
#endif
